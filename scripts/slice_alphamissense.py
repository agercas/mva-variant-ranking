#!/usr/bin/env python3
"""Slice the AlphaMissense bulk table down to a set of gene spans.

The full hg38 table is ~614 MB gzipped and ~71 M rows. Every analysis here wants
a few genes out of it, so this reduces it once to a few megabytes that the rest
of the pipeline can read instantly.

Args:
    --table: ``AlphaMissense_hg38.tsv.gz`` as downloaded.
    --genes: ``chrom<TAB>start<TAB>end<TAB>symbol``, no ``chr`` prefix (the
        Ensembl convention this repository uses). AlphaMissense is chr-prefixed;
        the mapping is applied here so no caller has to remember it.
    --out: the sliced table, each row prefixed with its gene symbol.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--genes", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    spans = []
    for line in args.genes.read_text().splitlines():
        chrom, start, end, gene = line.split("\t")
        spans.append((f"chr{chrom}", int(start), int(end), gene))

    kept = 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.table, "rt") as fh, args.out.open("w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _rest = line.split("\t", 2)
            position = int(pos)
            for span_chrom, start, end, gene in spans:
                if chrom == span_chrom and start <= position <= end:
                    out.write(f"{gene}\t{line}")
                    kept += 1
                    break
    print(f"{kept:,} rows -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
