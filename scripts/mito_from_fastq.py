#!/usr/bin/env python3
"""Call mitochondrial variants by streaming the raw FASTQs. Never stores them.

Why this exists
---------------
The supplied call set (the VCF) has no mitochondrial genome: zero MT records and
zero MT contig headers. That leaves a mitochondrial cause untested. The VCF does
not contain chrM, but the dataset does: the 85 GB of raw FASTQ holds the reads.

The obstacle is disk, not permission: 85 GB against ~21 GB free. So the reads are
**streamed and discarded**, never landed:

    HTTPS range-stream -> incremental gunzip -> k-mer prefilter -> minimap2 (mappy)
                       -> pileup accumulator -> a few hundred KB of counts

The k-mer prefilter is what makes this tractable. Aligning every read to MT would
be wasteful when >99.9% of them are nuclear; a 25-mer membership test against the
16,569 bp reference (both strands) rejects them for ~0.6 us each, and only
survivors reach the aligner.

🔴 **NUMTs are the standing caveat.** Nuclear mitochondrial insertions carry
MT-like sequence, so a bait pulls some of them too. Best practice is to align to
the whole genome and keep reads that map *best* to MT, which needs a reference
and RAM this box does not have. Mitigations, all reported rather than assumed:

* a high identity floor and a mapping-quality floor,
* coverage uniformity along the molecule (NUMT contamination is patchy, real MT
  is not),
* and haplogroup coherence -- a real mitochondrion yields a self-consistent set
  of known haplogroup-defining variants, whereas NUMT noise does not.

🔴 Reads only public references and the gated FASTQ. Nothing is transmitted:
the HTTPS traffic is *inbound* dataset download, exactly as ``hf download`` would
be, and no proband-derived coordinate is sent anywhere.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import zlib
from pathlib import Path

import mappy
import requests

REPO = Path(__file__).resolve().parent.parent
DATASET = "SageBio/mva-hackathon-2026-data"
BASE = f"https://huggingface.co/datasets/{DATASET}/resolve/main"

K = 25
MIN_KMER_HITS = 2
MIN_IDENTITY = 0.95
MIN_MAPQ = 20
MIN_ALN_LEN = 80

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def build_kmer_set(ref_seq: str, k: int = K) -> set[str]:
    """Every k-mer of the reference and its reverse complement."""
    kmers = {ref_seq[i : i + k] for i in range(len(ref_seq) - k + 1)}
    # The molecule is circular: bridge the origin so reads spanning it survive.
    wrapped = ref_seq[-(k - 1) :] + ref_seq[: k - 1]
    kmers |= {wrapped[i : i + k] for i in range(len(wrapped) - k + 1)}
    return kmers | {revcomp(x) for x in kmers}


def looks_mitochondrial(seq: str, kmers: set[str], k: int = K) -> bool:
    """Cheap membership test: does this read share enough k-mers with MT?"""
    hits = 0
    for start in range(0, len(seq) - k + 1, k):
        if seq[start : start + k] in kmers:
            hits += 1
            if hits >= MIN_KMER_HITS:
                return True
    return False


def stream_fastq_reads(url: str, token: str, byte_limit: int | None = None):
    """Yield read sequences from a remote gzipped FASTQ without storing it."""
    headers = {"Authorization": f"Bearer {token}"}
    if byte_limit:
        headers["Range"] = f"bytes=0-{byte_limit - 1}"
    # 🔴 These FASTQs are MULTI-MEMBER gzip (BGZF-style: many concatenated
    # streams). A single decompressobj stops at the end of the FIRST member and
    # then silently returns b"" forever -- which reads as "the file ended" and
    # cost a smoke test that found 177 reads in 400 MB. Roll a fresh
    # decompressor over each member, handing it the previous one's unused_data.
    decomp = zlib.decompressobj(zlib.MAX_WBITS | 16)
    pending = b""
    line_no = 0
    with requests.get(url, headers=headers, stream=True, timeout=600) as response:
        response.raise_for_status()
        for chunk in response.iter_content(1 << 22):
            while chunk:
                try:
                    text = decomp.decompress(chunk)
                except zlib.error:
                    return
                if text:
                    pending += text
                    *lines, pending = pending.split(b"\n")
                    for line in lines:
                        if line_no % 4 == 1:
                            yield line.decode("ascii", "replace")
                        line_no += 1
                if not decomp.eof:
                    break
                chunk = decomp.unused_data
                decomp = zlib.decompressobj(zlib.MAX_WBITS | 16)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mt-fasta", type=Path, required=True, help="rCRS / GRCh38 MT reference FASTA")
    ap.add_argument("--files", nargs="+", required=True, help="FASTQ filenames within the dataset repo")
    ap.add_argument(
        "--gb-per-file", type=float, default=0.0, help="stop after this many GB of each file (0 = all)"
    )
    ap.add_argument("--out", type=Path, default=REPO / "data/private/work/mito_pileup.json")
    args = ap.parse_args()

    token = os.environ.get("HF_TOKEN")
    if not token:
        print("HF_TOKEN not set")
        return 1

    ref_seq = "".join(
        line.strip() for line in args.mt_fasta.read_text().splitlines() if not line.startswith(">")
    )
    print(f"MT reference: {len(ref_seq):,} bp")
    kmers = build_kmer_set(ref_seq)
    print(f"k-mer bait  : {len(kmers):,} {K}-mers (both strands, origin-bridged)")

    aligner = mappy.Aligner(str(args.mt_fasta), preset="sr", n_threads=3)
    if not aligner:
        print("failed to build the minimap2 index")
        return 1

    pileup: dict[int, collections.Counter] = collections.defaultdict(collections.Counter)
    totals = {"reads": 0, "baited": 0, "aligned": 0, "bases": 0}

    byte_limit = int(args.gb_per_file * 1e9) if args.gb_per_file else None
    for name in args.files:
        print(f"\n--- streaming {name} ---", flush=True)
        for seq in stream_fastq_reads(f"{BASE}/{name}", token, byte_limit):
            totals["reads"] += 1
            if totals["reads"] % 20_000_000 == 0:
                print(
                    f"    {totals['reads']:,} reads | baited {totals['baited']:,} | "
                    f"aligned {totals['aligned']:,}",
                    flush=True,
                )
            if not looks_mitochondrial(seq, kmers):
                continue
            totals["baited"] += 1
            for hit in aligner.map(seq):
                if hit.mapq < MIN_MAPQ or (hit.q_en - hit.q_st) < MIN_ALN_LEN:
                    continue
                if hit.mlen / max(hit.blen, 1) < MIN_IDENTITY:
                    continue
                totals["aligned"] += 1
                query = seq if hit.strand == 1 else revcomp(seq)
                qpos, rpos = hit.q_st, hit.r_st
                for length, op in hit.cigar:  # op: 0=M 1=I 2=D 3=N 4=S
                    if op == 0:
                        for i in range(length):
                            pileup[rpos + i][query[qpos + i]] += 1
                        qpos += length
                        rpos += length
                        totals["bases"] += length
                    elif op in (1, 4):
                        qpos += length
                    elif op in (2, 3):
                        for i in range(length):
                            pileup[rpos + i]["*"] += 1
                        rpos += length
                break  # primary hit only

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(
            {
                "totals": totals,
                "reference_length": len(ref_seq),
                "pileup": {str(k): dict(v) for k, v in sorted(pileup.items())},
            }
        )
    )
    covered = len(pileup)
    depth = totals["bases"] / max(covered, 1)
    print(f"\nreads streamed : {totals['reads']:,}")
    print(f"k-mer baited   : {totals['baited']:,}  ({100 * totals['baited'] / max(totals['reads'], 1):.4f}%)")
    print(f"aligned to MT  : {totals['aligned']:,}")
    print(f"positions covered: {covered:,} / {len(ref_seq):,}  ({100 * covered / len(ref_seq):.1f}%)")
    print(f"mean depth     : {depth:.0f}x")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
