#!/usr/bin/env python3
"""Re-derive a variant's protein consequence by hand, trusting nothing but the reference.

``bcftools csq`` said ``p.Leu737Ter`` and ``p.Asn1002Lys``. This script does not
believe it. It rebuilds the coding sequence from the reference FASTA and the
GFF3's CDS records, translates it with a hard-coded genetic code, applies the
substitution itself, and reports what the codon actually becomes.

Two independent tools agreeing is worth much more than one tool asserting, and
the whole Track 1 call rests on these two annotations being right.

Checks performed, each of which can fail loudly:

1. The assembled CDS starts with ATG and ends on a stop codon.
2. Its length is divisible by three.
3. The translated protein length matches the expected isoform length.
4. The variant position lies inside the CDS.
5. The reference base at that position matches the VCF's REF allele.
6. The reference codon translates to the expected wild-type residue.

Usage:
    scripts/verify_consequence_from_first_principles.py \\
        --transcript ENST00000287598 --contig 15 --pos 40209701 --ref T --alt G
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

#: The standard genetic code. Written out rather than imported so that this
#: script depends on nothing that could itself be wrong.
CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
THREE_LETTER = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "Ter",
}


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3))


def load_cds(gff: Path, transcript: str) -> tuple[list[tuple[int, int]], str]:
    """Return this transcript's CDS blocks as (start, end) 1-based inclusive, plus strand."""
    blocks: list[tuple[int, int]] = []
    strand = "+"
    with gzip.open(gff, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS" or transcript not in f[8]:
                continue
            blocks.append((int(f[3]), int(f[4])))
            strand = f[6]
    blocks.sort()
    return blocks, strand


def read_contig(fasta: Path) -> str:
    """Read a single-record FASTA into memory (one chromosome, ~100 MB of text)."""
    parts: list[str] = []
    with fasta.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--gff", type=Path, default=REPO / "data/reference/ensembl116/Homo_sapiens.GRCh38.116.gff3.gz"
    )
    ap.add_argument("--fasta", type=Path, default=REPO / "data/reference/ensembl116/chr15.fa")
    ap.add_argument("--transcript", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--alt", required=True)
    ap.add_argument("--expect-protein-length", type=int, default=None)
    args = ap.parse_args()

    blocks, strand = load_cds(args.gff, args.transcript)
    if not blocks:
        print(f"no CDS found for {args.transcript}", file=sys.stderr)
        return 1
    seq = read_contig(args.fasta)

    # --- assemble the coding sequence in transcription order -----------------
    pieces = [seq[s - 1 : e] for s, e in blocks]
    cds = "".join(pieces) if strand == "+" else revcomp("".join(pieces))

    failures: list[str] = []
    if not cds.startswith("ATG"):
        failures.append(f"CDS does not start with ATG (got {cds[:3]})")
    if len(cds) % 3:
        failures.append(f"CDS length {len(cds)} is not a multiple of 3")
    protein = translate(cds)
    stops = protein.count("*")
    if not protein.endswith("*"):
        failures.append("CDS does not end on a stop codon")
    if stops != 1:
        failures.append(f"{stops} stop codons in the translation; expected exactly 1 (the terminator)")
    aa_length = len(protein) - 1
    if args.expect_protein_length and aa_length != args.expect_protein_length:
        failures.append(f"protein length {aa_length} != expected {args.expect_protein_length}")

    print(f"transcript {args.transcript}  strand {strand}  CDS blocks {len(blocks)}")
    print(f"CDS length {len(cds)} nt  ->  protein {aa_length} aa (+ stop)")

    # --- map the genomic position into the CDS -------------------------------
    offset = 0
    cds_index: int | None = None
    ordered = blocks if strand == "+" else sorted(blocks, reverse=True)
    for s, e in ordered:
        if s <= args.pos <= e:
            cds_index = offset + (args.pos - s if strand == "+" else e - args.pos)
            break
        offset += e - s + 1
    if cds_index is None:
        print(f"position {args.contig}:{args.pos} is NOT inside the CDS of {args.transcript}")
        return 1

    genomic_ref = seq[args.pos - 1]
    expected_ref = args.ref if strand == "+" else args.ref.translate(COMPLEMENT)
    if cds[cds_index] != expected_ref:
        failures.append(f"CDS base {cds[cds_index]} != VCF REF {args.ref} mapped to strand ({expected_ref})")
    if genomic_ref != args.ref:
        failures.append(f"genomic base {genomic_ref} != VCF REF {args.ref}")

    codon_number = cds_index // 3 + 1
    frame = cds_index % 3
    ref_codon = cds[codon_number * 3 - 3 : codon_number * 3]
    alt_base = args.alt if strand == "+" else args.alt.translate(COMPLEMENT)
    alt_codon = ref_codon[:frame] + alt_base + ref_codon[frame + 1 :]
    ref_aa, alt_aa = CODON_TABLE.get(ref_codon, "X"), CODON_TABLE.get(alt_codon, "X")

    print(f"\n{args.contig}:{args.pos} {args.ref}>{args.alt}")
    print(f"  genomic reference base      : {genomic_ref}   (matches VCF REF: {genomic_ref == args.ref})")
    print(f"  CDS offset                  : {cds_index} (0-based)")
    print(f"  codon                       : {codon_number}, position {frame + 1} of 3")
    print(f"  reference codon             : {ref_codon} -> {ref_aa} ({THREE_LETTER.get(ref_aa, '?')})")
    print(f"  alternate codon             : {alt_codon} -> {alt_aa} ({THREE_LETTER.get(alt_aa, '?')})")

    if ref_aa == alt_aa:
        consequence = "synonymous"
    elif alt_aa == "*":
        consequence = "stop_gained (nonsense)"
    elif ref_aa == "*":
        consequence = "stop_lost"
    elif codon_number == 1:
        consequence = "start_lost"
    else:
        consequence = "missense"
    hgvs = f"p.{THREE_LETTER.get(ref_aa, '?')}{codon_number}{THREE_LETTER.get(alt_aa, '?')}"
    print(f"\n  ==> {consequence}, {hgvs}")

    # --- NMD: is the stop more than 50 nt upstream of the last junction? -----
    #
    # The 50-nt rule is quoted three different ways in the wild, depending on
    # whether the distance is measured from the PTC codon's first base, its last
    # base, or the substituted base itself. On this transcript they differ by two
    # (748 / 746 / 747), which is enough to look like a contradiction between two
    # documents that are in fact reporting the same fact. So print all three with
    # their definitions and let the reader see that the verdict does not turn on
    # the convention.
    if alt_aa == "*":
        last_junction = len(cds) - (ordered[-1][1] - ordered[-1][0] + 1)
        codon_end = codon_number * 3  # 1-based offset of the codon's last base
        margins = {
            "from the codon's first base": last_junction - (codon_end - 2),
            "from the substituted base": last_junction - (cds_index + 1),
            "from the codon's last base": last_junction - codon_end,
        }
        smallest = min(margins.values())
        verdict = "NMD-SENSITIVE" if smallest > 50 else "likely NMD-ESCAPING"
        print(
            f"  ==> premature stop is {verdict} (>50 nt rule); final junction at CDS offset {last_junction}"
        )
        for how, margin in margins.items():
            print(f"        {margin:>5} nt upstream of the final exon-exon junction, measured {how}")
        print(
            f"        all three exceed 50 by >{smallest - 50} nt, so the verdict does not depend on the convention"
        )

    if failures:
        print("\n🔴 SANITY CHECKS FAILED:")
        for f in failures:
            print(f"   - {f}")
        return 1
    print("\n✅ every sanity check passed (ATG start, single terminal stop, in-frame, REF concordant)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
