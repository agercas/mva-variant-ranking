#!/usr/bin/env python3
"""Build and validate the Track 1 predictions CSV.

Usage::

    python3 scripts/build_track1_submission.py --out <path>.csv

``--out`` is required, so a bare run cannot overwrite a file by accident.

Row 1 is the compound-heterozygous pair, in ONE row. Two rows would be scored as
a partial match and capped at 50 rank points. Rows 2-4 sit strictly below it: a
full match at rank 1 suppresses partial credit and F-max is reached at the top
threshold, so the lower rows cost nothing if row 1 is right and keep 50 points
if it is only half right.

This script writes a file. It does not upload anything. Spending one of the six
submission slots is a human decision.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mva.genomics.coordinates import Variant
from mva.genomics.submission import Prediction, validate, write
from mva.scoring import score

REPO = Path(__file__).resolve().parent.parent

#: p.Leu737Ter on ENST00000287598 / NM_001211.6. Exon 17 of 23, so the premature
#: stop sits far upstream of the last exon-exon junction and is expected to trigger
#: nonsense-mediated decay; any protein that escaped would still lose the whole
#: kinase domain (UniProt O60566, residues 766-1050). ClinVar 533901
#: Pathogenic/Likely_pathogenic for mosaic variegated aneuploidy syndrome 1.
#: gnomAD v4.1 exomes AF 7.87e-05, 0 homozygotes.
NONSENSE = Variant("15", 40209701, "T", "G")
#: p.Asn1002Lys, last exon, inside the kinase domain. Not in ClinVar under this
#: nucleotide change (ClinVar 4600147 is the T>A change that gives the same protein
#: change, classified VUS). A singleton in gnomAD v4.1 exomes: AC=1 of AN=1,461,878.
MISSENSE = Variant("15", 40220612, "T", "G")
#: Deep intronic, 1,026 bp from the nearest exon boundary, not protein-altering.
#: Kept only as a hedge in case the second hit is not the missense.
INTRONIC = Variant("15", 40216470, "A", "G")
#: Secondary finding, unrelated to the primary diagnosis. LZTR1 p.Tyr748Ter
#: (ENST00000646124), stop-gain expected to trigger nonsense-mediated decay,
#: ClinVar 1409252 Pathogenic/Likely_pathogenic.
LZTR1_SECONDARY = Variant("22", 20996720, "C", "G")


def predictions() -> list[Prediction]:
    """The five rows, best first.

    Row 1 carries almost all of the belief. Across a 170 kb window that covers the
    upstream region, the whole BUB1B gene and the BUB1B-PAK6 readthrough tail,
    these are the only two protein-altering variants the child has, so there is
    exactly one coding compound-heterozygous pair available.

    ``epcr`` is meant to be a probability. We cannot calculate one, so the values only
    sort the rows, best first. The organizers' scorer uses them only to sort rows and
    to set the F-max thresholds, so their size does not change the score.

    Rows 2 and 3 are ordered by how plausible we find each alternative. Row 2
    (the stop-gain alone) covers a second hit that this call set cannot show. Row 3
    (stop-gain plus the deep intronic variant) is less likely, because SpliceAI
    predicts no splicing effect for that variant. Ranks 2 and 3 pay the same rank
    points, so the order only reflects belief.
    """
    return [
        Prediction(
            "PROBAND01",
            NONSENSE,
            0.9,
            "primary",
            MISSENSE,
            "BUB1B compound heterozygote (NM_001211.6 / ENST00000287598). Allele 1: p.Leu737Ter, "
            "stop-gain in exon 17 of 23, expected to trigger nonsense-mediated decay; ClinVar 533901 "
            "Pathogenic/Likely pathogenic for mosaic variegated aneuploidy syndrome 1; gnomAD AF 7.9e-05. "
            "Allele 2: p.Asn1002Lys, missense in the kinase domain, seen once in 1.46 million gnomAD "
            "alleles; not in ClinVar under this nucleotide change. These are the only two "
            "protein-altering variants in the gene. Phase cannot be determined from this data.",
        ),
        Prediction(
            "PROBAND01",
            NONSENSE,
            0.6,
            "primary",
            notes="Stop-gain allele alone, in case the second hit is something this call set cannot "
            "show (for example a balanced rearrangement). Exon-level deletions are excluded: all 23 "
            "coding exons have 0.77 to 1.44 times the median read depth.",
        ),
        Prediction(
            "PROBAND01",
            NONSENSE,
            0.4,
            "primary",
            INTRONIC,
            "Stop-gain allele plus the only other novel (absent from dbSNP) heterozygous variant in BUB1B, a deep intronic "
            "change 1,026 bp from the nearest exon. Less likely than row 2: SpliceAI predicts no "
            "splicing effect (0.00).",
        ),
        Prediction("PROBAND01", MISSENSE, 0.2, "primary", notes="Missense allele alone."),
        Prediction(
            "PROBAND01",
            LZTR1_SECONDARY,
            0.1,
            "secondary",
            notes="SECONDARY FINDING, unrelated to the primary diagnosis. LZTR1 p.Tyr748Ter "
            "(ENST00000646124): stop-gain expected to trigger nonsense-mediated decay; ClinVar 1409252 "
            "Pathogenic/Likely pathogenic. Heterozygous loss of LZTR1 function is associated with "
            "schwannomatosis (autosomal dominant, incomplete penetrance), for which clinical follow-up "
            "can be considered. It does not explain this condition: Noonan syndrome 10 is caused by "
            "dominant missense variants and Noonan syndrome 2 needs variants in both copies of LZTR1. "
            "Needs confirmation in an accredited laboratory and review by the treating team.",
        ),
    ]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True, help="where to write the CSV")
    args = ap.parse_args()

    preds = predictions()
    report = validate(preds)
    print(report.render())
    if not report.ok:
        return 1

    path = write(preds, args.out)
    # --out may point outside the repo (a temporary directory in a test, a scratch
    # path for a comparison), and relative_to() raises for those.
    try:
        shown: Path | str = path.relative_to(REPO)
    except ValueError:
        shown = path
    print(f"\nwrote {shown}\n")

    # Offline dry run against a copy of the organizers' scoring rules, under the
    # three answers we can name. This shows the ladder behaves as intended. It
    # proves nothing about which answer is true.
    scenarios = {
        "our pair is the answer": frozenset({NONSENSE, MISSENSE}),
        "stop-gain + intronic": frozenset({NONSENSE, INTRONIC}),
        "stop-gain alone": frozenset({NONSENSE}),
    }
    print(f"{'scenario':<26} {'rank points':>12} {'F-max':>8}")
    print("-" * 48)
    for label, truth in scenarios.items():
        result = score(preds, truth)
        print(f"{label:<26} {result.rank_points:>12.1f} {result.f_max:>8.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
