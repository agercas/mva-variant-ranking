#!/usr/bin/env python3
"""Where do the BUB1B missense alleles fall in AlphaMissense's own scores?

AlphaMissense (Cheng et al., Science 2023) has already scored every possible
amino-acid change in every human protein, and the tables are public. So this
script does not run a model. It reads the table for BUB1B (UniProt O60566) and
answers four questions about the child's missense allele, p.Asn1002Lys:

1. What score does it have, and how does that rank among all BUB1B changes?
2. Is that score special for lysine, or is every change at residue 1002 scored high?
3. Which ClinGen evidence strength would the score map to if the tool had been
   chosen in advance? (Bands: Bergquist et al., Genet Med 2025, Table 1.)
4. What does the tool say about BUB1B missense alleles whose effect on BubR1
   protein was measured in the laboratory (Suijkerbuijk et al., Cancer Res 2010)?

No data from the child is read. The only input is the public table.

Usage::

    python3 scripts/alphamissense_bub1b.py --download data/reference/alphamissense
    python3 scripts/alphamissense_bub1b.py --dir data/reference/alphamissense

``--download`` fetches the two per-protein files from AlphaFold DB. The tables
are not stored in this repository.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import statistics
import urllib.request
from collections import defaultdict
from pathlib import Path

BASE = "https://alphafold.ebi.ac.uk/files/"
FILES = ("AF-O60566-F1-aa-substitutions.csv", "AF-O60566-F1-hg38.csv")

#: The AlphaMissense developers' own boundary between "ambiguous" and "likely pathogenic".
DEVELOPER_PATHOGENIC_CUTOFF = 0.564

#: ClinGen-calibrated evidence bands for AlphaMissense, Bergquist et al. 2025
#: (PMID 40084623, PMC12208618), Table 1. (lower, upper, label, ACMG points).
#: Scores between 0.170 and 0.791 are indeterminate (no criterion).
AM_BANDS: list[tuple[float, float, str, int]] = [
    (0.000, 0.070, "BP4 (strength -3)", -3),
    (0.071, 0.099, "BP4_moderate", -2),
    (0.100, 0.169, "BP4_supporting", -1),
    (0.170, 0.791, "indeterminate", 0),
    (0.792, 0.905, "PP3_supporting", 1),
    (0.906, 0.971, "PP3_moderate", 2),
    (0.972, 0.989, "PP3 (strength +3)", 3),
    (0.990, 1.000, "PP3_strong", 4),
]

QUERY = "N1002K"

#: Reference missense alleles with a laboratory measurement.
#: Source: Suijkerbuijk et al., Cancer Res 2010 (PMID 20516114). "low" means the BubR1
#: protein was less abundant than wild type (patient cells or transfected cells);
#: "normal" means indistinguishable from wild type.
REFERENCE_ALLELES: list[tuple[str, bool, str]] = [
    ("R727C", True, "low abundance in patient cells"),
    ("R814H", True, "low abundance in patient cells"),
    ("L844F", True, "poorly expressed, not functional"),
    ("I909T", True, "poorly expressed; function restored when overexpressed"),
    ("L1012P", True, "poorly expressed; function restored when overexpressed"),
    ("Q921H", False, "indistinguishable from wild type"),
]


def band(score: float) -> tuple[str, int]:
    for lo, hi, label, points in AM_BANDS:
        if lo <= score <= hi:
            return label, points
    return "indeterminate", 0  # gaps between printed bands (rounding) carry no criterion


def share_below(sorted_scores: list[float], value: float) -> float:
    return 100.0 * bisect.bisect_left(sorted_scores, value) / len(sorted_scores)


def residue(variant: str) -> int:
    return int(variant[1:-1])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, default=Path("data/reference/alphamissense"))
    ap.add_argument("--download", type=Path, help="fetch the two public tables into this directory first")
    args = ap.parse_args()

    if args.download:
        args.download.mkdir(parents=True, exist_ok=True)
        for name in FILES:
            req = urllib.request.Request(BASE + name, headers={"User-Agent": "Mozilla/5.0 research-script"})
            (args.download / name).write_bytes(urllib.request.urlopen(req, timeout=120).read())
        args.dir = args.download

    aa_rows = list(csv.DictReader((args.dir / FILES[0]).open()))
    nt_rows = list(csv.DictReader((args.dir / FILES[1]).open()))
    score = {r["protein_variant"]: float(r["am_pathogenicity"]) for r in aa_rows}
    q = score[QUERY]

    print("AlphaMissense scores for BUB1B (UniProt O60566). Public table, no data from the child.")
    print(f"amino-acid changes in the table: {len(score):,}  (1,050 residues x 19)")
    print(f"changes reachable by one nucleotide change: {len(nt_rows):,} nucleotide-level rows\n")

    # 1. The score and where it ranks
    all_sorted = sorted(score.values())
    nt_sorted = sorted(float(r["am_pathogenicity"]) for r in nt_rows)
    kin_nt = sorted(
        float(r["am_pathogenicity"]) for r in nt_rows if 766 <= residue(r["protein_variant"]) <= 1050
    )
    chrom_pos = [
        r for r in nt_rows if r["protein_variant"] == QUERY and r["POS"] == "40220612" and r["ALT"] == "G"
    ]
    print("1. THE SCORE")
    print(
        f"   p.Asn1002Lys = {q}  (class: {next(r['am_class'] for r in aa_rows if r['protein_variant'] == QUERY)})"
    )
    print(
        f"   the nucleotide-level table lists chr15:40220612 T>G as {QUERY}: "
        f"{'yes, score ' + chrom_pos[0]['am_pathogenicity'] if chrom_pos else 'NOT FOUND'}"
    )
    print(
        f"   {share_below(nt_sorted, q):.1f}th percentile of the {len(nt_sorted):,} one-nucleotide changes in BUB1B"
    )
    print(
        f"   {share_below(kin_nt, q):.1f}th percentile of the {len(kin_nt):,} one-nucleotide changes in residues 766-1050"
    )
    higher = len(all_sorted) - bisect.bisect_right(all_sorted, q)
    print(
        f"   {share_below(all_sorted, q):.1f}th percentile of all {len(all_sorted):,} changes "
        f"({higher:,} score higher)\n"
    )

    # 2. Is the score about lysine, or about the position?
    site = sorted(((s, v[-1]) for v, s in score.items() if residue(v) == 1002), reverse=True)
    print("2. THE SAME-SITE LADDER (all 19 changes at residue 1002, highest first)")
    for rank, (s, aa) in enumerate(site, 1):
        print(f"   {rank:>2}. {aa}  {s:.4f}{'   <-- lysine' if aa == 'K' else ''}")
    above = sum(s >= DEVELOPER_PATHOGENIC_CUTOFF for s, _ in site)
    print(
        f"   {above} of 19 changes at residue 1002 score above the developers' cutoff of "
        f"{DEVELOPER_PATHOGENIC_CUTOFF}; median {statistics.median(s for s, _ in site):.3f}; "
        f"lysine is rank {[a for _, a in site].index('K') + 1} of 19"
    )
    by_site: dict[int, list[float]] = defaultdict(list)
    for v, s in score.items():
        by_site[residue(v)].append(s)
    means = {k: statistics.mean(v) for k, v in by_site.items()}
    print(
        f"   mean over the 19 changes = {means[1002]:.2f}, higher than "
        f"{100.0 * sum(m < means[1002] for m in means.values()) / len(means):.0f}% of residues in BUB1B"
    )
    print(
        "   Reading: the position is scored as intolerant to almost any change. Most of the "
        "0.92 comes from the position, not from lysine.\n"
    )

    # 3. ClinGen-calibrated band
    label, points = band(q)
    print("3. WHAT THE CLINGEN-CALIBRATED BANDS WOULD SAY (Bergquist 2025, Table 1)")
    print(f"   {q} falls in the band '{label}' ({points:+d} points).")
    print(
        "   The developers' own cutoff (0.564) is below the lowest evidence band (0.792), "
        "so the tool's 'likely pathogenic' label is not an ACMG strength."
    )
    print(
        "   We do NOT use this as a criterion: ClinGen asks for one tool to be chosen before the "
        "variant is looked at, and we did not do that. See the report, section 5.\n"
    )

    # 4. The tool against alleles with a measured effect
    print("4. THE TOOL ON ALLELES WITH A LABORATORY MEASUREMENT (Suijkerbuijk 2010, PMID 20516114)")
    print(f"   {'allele':<8}{'score':>8}  {'class':<8}{'ClinGen band':<18}{'laboratory result'}")
    cls = {r["protein_variant"]: r["am_class"] for r in aa_rows}
    low = [v for v, is_low, _ in REFERENCE_ALLELES if is_low]
    missed = [v for v in low if score[v] < DEVELOPER_PATHOGENIC_CUTOFF]
    for v, _, note in REFERENCE_ALLELES:
        s = score[v]
        lab, _ = band(s)
        print(f"   {v:<8}{s:>8.4f}  {cls[v]:<8}{lab:<18}{note}")
    print(
        f"   {len(missed)} of the {len(low)} alleles with low protein abundance score below the developers' "
        f"cutoff ({', '.join(missed)})."
    )
    print(
        "   Reading: a high score supports damage, but a low score would not have ruled it out. "
        "The tool does not say how much BubR1 protein is left."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
