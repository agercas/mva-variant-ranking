#!/usr/bin/env python3
"""Classify the two BUB1B alleles with the ACMG/AMP criteria, on the point scale.

A diagnostic laboratory would put each allele through the ACMG/AMP criteria
(Richards et al. 2015). This script does that and prints, for each criterion,
whether we counted it and why. It also answers two questions a reader will ask:

* Would a different missense predictor change the class of allele 2?
* What further evidence would be needed to lift allele 2 out of "uncertain"?

Points (Tavtigian et al., Hum Mutat 2020): Supporting 1, Moderate 2, Strong 4,
Very Strong 8. Pathogenic >= 10, Likely pathogenic 6-9, Uncertain significance
(VUS) 0-5, Likely benign -1 to -6, Benign <= -7.

Every number printed here is computed from the lists below. The evidence behind
each criterion is in ``results/`` (see the README of the submission repository).
"""

from __future__ import annotations

from dataclasses import dataclass

WEIGHT = {"supporting": 1, "moderate": 2, "strong": 4, "very_strong": 8}


@dataclass(frozen=True)
class Criterion:
    code: str
    strength: str
    applied: bool
    reason: str

    @property
    def points(self) -> int:
        return WEIGHT[self.strength] if self.applied else 0


def classify(points: int) -> str:
    if points >= 10:
        return "PATHOGENIC"
    if points >= 6:
        return "LIKELY PATHOGENIC"
    if points >= 0:
        return "UNCERTAIN SIGNIFICANCE (VUS)"
    if points >= -6:
        return "LIKELY BENIGN"
    return "BENIGN"


ALLELE1 = [
    Criterion(
        "PVS1",
        "very_strong",
        True,
        "Stop-gain p.Leu737Ter in exon 17 of 23. The stop is 747 nucleotides "
        "upstream of the last exon-exon junction (the usual rule is more than 50), "
        "so nonsense-mediated decay is expected. Loss of function of one copy is the "
        "mechanism seen in patients: the first patients described carried a truncating "
        "variant on one copy (Hanks 2004), and most pathogenic BUB1B variants in ClinVar "
        "are nonsense, frameshift or splice-site changes. The position was checked against "
        "the Ensembl transcript.",
    ),
    Criterion(
        "PM2",
        "supporting",
        False,
        "PM2 asks for absence or extreme rarity. This allele is rare but "
        "seen 120 times in gnomAD v4.1 (exomes and genomes), never homozygous, so we "
        "do not count it. Counting it would give 9 points; the class is the same.",
    ),
    Criterion(
        "PM3",
        "moderate",
        False,
        "It needs the variant in trans with a pathogenic allele, and "
        "phase cannot be measured in this data: the two variants are 10,911 bp "
        "apart and the longest DNA fragment seen at this locus is 800 bp.",
    ),
    Criterion("PS3", "strong", False, "No laboratory test of this allele has been published."),
    Criterion(
        "PP5",
        "supporting",
        False,
        "ClinVar lists it as Pathogenic/Likely pathogenic with several "
        "submitters and no conflicts, which supports our reading, but ClinGen no "
        "longer recommends PP5 as a criterion.",
    ),
]

ALLELE2 = [
    Criterion(
        "PM2",
        "supporting",
        True,
        "Seen once in gnomAD v4.1 exomes (1 of 1,461,878 alleles), never homozygous.",
    ),
    Criterion(
        "PP3",
        "supporting",
        False,
        "ClinGen asks for ONE predictor, chosen before the variant is "
        "looked at. We did not choose one in advance. The predictors also disagree "
        "(see the table below). Counting any one of them would not change the "
        "class: the table shows every choice leaves the allele uncertain.",
    ),
    Criterion(
        "PM1",
        "moderate",
        False,
        "PM1 needs a mutational hot spot or a site known to be critical for function. "
        "N1002 is in the kinase domain (UniProt 766-1050), which BubR1 uses as a "
        "scaffold (it is generally regarded as a pseudokinase), and no laboratory result "
        "marks this residue as critical. In the AlphaFold model it is mostly buried "
        "(relative solvent accessibility 0.14) and conserved (asparagine in 111 of 115 "
        "aligned vertebrate sequences). That makes a structural role plausible. It is "
        "not the experimental evidence PM1 asks for.",
    ),
    Criterion("PM3", "moderate", False, "Same reason as allele 1: phase cannot be measured here."),
    Criterion("PS3", "strong", False, "No laboratory test of p.Asn1002Lys exists."),
    Criterion(
        "PP4",
        "supporting",
        False,
        "PP4 needs a phenotype specific to a gene with a single "
        "genetic cause. Mosaic variegated aneuploidy has several known genes.",
    ),
]

#: What each calibrated predictor says about p.Asn1002Lys. The scores come from
#: results/insilico-panel.txt and results/alphamissense-bub1b.txt. The bands are the
#: ClinGen-calibrated ones: Pejaver et al. 2022 (PMID 36413997) for REVEL, VEST4,
#: PolyPhen-2 and SIFT; Bergquist et al. 2025 (PMID 40084623) for AlphaMissense.
#: (tool, score, criterion it would give, ACMG points it would add)
#: The VEST4 score is the value the UCSC track gave during the analysis. It does not reproduce
#: against dbNSFP 4.8a (0.885 on the MANE transcript, which is PP3_moderate); see
#: results/RECORDED-CHECKS.md. The class of allele 2 does not depend on it.
PREDICTORS: list[tuple[str, str, str, int]] = [
    ("REVEL", "0.472", "indeterminate (no criterion)", 0),
    ("VEST4", "0.433", "BP4_supporting", -1),
    ("PolyPhen-2", "0.997", "PP3_supporting", 1),
    ("SIFT", "0.000", "PP3_moderate", 2),
    ("AlphaMissense", "0.9229", "PP3_moderate", 2),
]


def report(name: str, desc: str, criteria: list[Criterion]) -> int:
    print(f"\n{'=' * 74}\n{name}\n  {desc}\n{'=' * 74}")
    total = sum(c.points for c in criteria)
    for c in criteria:
        mark = f"COUNTED  +{c.points}" if c.applied else "not counted"
        print(f"\n  {c.code:<5} [{c.strength}] {mark}")
        words, line = c.reason.split(), "        "
        for w in words:
            if len(line) + len(w) > 76:
                print(line.rstrip())
                line = "        "
            line += w + " "
        print(line.rstrip())
    print(f"\n  {'-' * 70}\n  TOTAL {total} points  ->  {classify(total)}")
    return total


t1 = report("ALLELE 1 - chr15:40209701 T>G", "c.2210T>G  p.Leu737Ter  (stop-gain)", ALLELE1)
t2 = report("ALLELE 2 - chr15:40220612 T>G", "c.3006T>G  p.Asn1002Lys  (missense)", ALLELE2)

print(f"\n{'=' * 74}\nWOULD A DIFFERENT PREDICTOR CHANGE THE CLASS OF ALLELE 2?\n{'=' * 74}")
print("  Total points for allele 2 if each tool had been the single tool chosen in advance:\n")
print(f"  {'tool':<15}{'score':>8}  {'criterion':<30}{'total':>6}  class")
for tool, score, crit, pts in PREDICTORS:
    total = t2 + pts
    print(f"  {tool:<15}{score:>8}  {crit:<30}{total:>6}  {classify(total)}")
print(
    "\n  Whatever the choice, allele 2 stays a VUS (highest total: "
    f"{max(t2 + p for *_, p in PREDICTORS)} points; Likely pathogenic starts at 6)."
)

print(f"\n{'=' * 74}\nWHAT WOULD LIFT ALLELE 2 OUT OF 'UNCERTAIN'\n{'=' * 74}")
print(f"  Allele 2 starts at {t2} point(s). Each line adds to that base.")
print("  PM3 uses the 2015 default weights; ClinGen's newer per-case PM3 scoring is not applied.\n")

#: (label, points added). Strengths are the defaults unless the label says otherwise.
#: A functional assay reaches Very Strong only if it is validated to a standard
#: (Brnich et al. 2019: a demonstrated OddsPath), which is a high bar.
SINGLE = [
    ("a calibrated predictor, chosen in advance, at Supporting (PP3)", WEIGHT["supporting"]),
    ("a calibrated predictor, chosen in advance, at Moderate (PP3)", WEIGHT["moderate"]),
    ("a calibrated predictor, chosen in advance, at Strong (PP3)", WEIGHT["strong"]),
    ("phase shown to be in trans (PM3 Moderate)", WEIGHT["moderate"]),
    ("phase in trans, plus a second unrelated patient (PM3 Strong)", WEIGHT["strong"]),
    ("a laboratory test showing less BubR1 protein (PS3 Strong)", WEIGHT["strong"]),
    ("a laboratory test validated to Very Strong (PS3 Very Strong)", WEIGHT["very_strong"]),
]
COMBINED = [
    ("phase (PM3 Strong) + a predictor at Supporting", WEIGHT["strong"] + WEIGHT["supporting"]),
    ("phase (PM3 Strong) + a laboratory test (PS3 Strong)", WEIGHT["strong"] + WEIGHT["strong"]),
    ("a predictor at Strong + a laboratory test (PS3 Strong)", WEIGHT["strong"] + WEIGHT["strong"]),
]

print("  ONE further piece of evidence:")
for label, pts in SINGLE:
    new = t2 + pts
    print(f"    +{pts}  {label:<64} -> {new:>2} points  {classify(new)}")

print("\n  TWO further pieces of evidence:")
for label, pts in COMBINED:
    new = t2 + pts
    print(f"    +{pts}  {label:<64} -> {new:>2} points  {classify(new)}")

default_single = [lbl for lbl, pts in SINGLE if t2 + pts >= 6 and "Very Strong" not in lbl]
print(f"\n  {'-' * 70}\n  Result:")
if not default_single:
    print("    No single piece of evidence at its default strength is enough. Phase alone")
    print("    is not enough, and neither is a laboratory test. The routes that work are a")
    print("    laboratory test validated to Very Strong, or two pieces together, for")
    print("    example phase at Strong with a predictor or a laboratory test (see above).")
else:
    print(f"    Single pieces that are enough at default strength: {', '.join(default_single)}")
print("    Phase is still the most informative single experiment, and it needs parents'")
print("    samples or long reads. This data set cannot supply it (longest fragment seen")
print("    at the locus 800 bp; the variants are 10,911 bp apart).")
