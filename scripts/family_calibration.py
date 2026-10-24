#!/usr/bin/env python3
"""Calibrate a missense predictor on the spindle-assembly-checkpoint gene family.

Why this exists
---------------
The ClinGen-calibrated predictors give no criterion for *BUB1B*'s best-known
MVA1 missense allele: REVEL scores `p.Leu1012Pro` at 0.500, which earns no
criterion in either direction. No Variant Curation Expert Panel exists for this
gene, so the generic genome-wide thresholds are all any laboratory has.

The generic thresholds come from pooling every gene in the genome. This script
asks the obvious next question: **do they look different if you pool the gene
family instead?** Spindle-checkpoint genes share a mechanism (mitotic checkpoint
complex function, dosage-sensitive in the same way), which is a tighter and more
defensible pooling assumption than "all genes".

🔴 **This cannot license PP3, and it is not used to.** ClinGen requires the tool
and threshold to be **pre-specified**, and this calibration was built *after*
seeing our variant's score. It is reported as a methods contribution and as a
measurement of how badly the generic thresholds fit this family -- never as
evidence for our own call. The report claims no computational criterion.

Method
------
Positive set: ClinVar Pathogenic/Likely-pathogenic missense across the family.
Negative set: ClinVar Benign/Likely-benign missense across the family.
Neither is chosen by us.

For a threshold *t*, the positive likelihood ratio is

    LR+(t) = P(score >= t | pathogenic) / P(score >= t | benign)

mapped to evidence strength by the Tavtigian OddsPath thresholds (prior 0.10):
supporting 2.08, moderate 4.33, strong 18.7, very strong 350.

With a small positive set, a point estimate is incomplete on its own, so the
script reports a **bootstrap interval** over the positives and, where no benign
variant clears the threshold, uses the **Rule of Three** (0 events in *n* trials
implies a rate below 3/*n* with ~95% confidence) rather than dividing by zero.
"""

from __future__ import annotations

import argparse
import random
import statistics
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CLINVAR = REPO / "data/reference/clinvar/clinvar.vcf.gz"

#: Tavtigian et al. 2018/2020 OddsPath thresholds at a prior of 0.10.
ODDS_PATH = ((350.0, "very strong"), (18.7, "strong"), (4.33, "moderate"), (2.08, "supporting"))


def strength_for(lr: float) -> str:
    for bound, label in ODDS_PATH:
        if lr >= bound:
            return f"PP3_{label}"
    return "below supporting"


def load_scores(path: Path) -> dict[tuple[str, int, str, str], tuple[str, float]]:
    out = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        gene, chrom, pos, ref, alt, score = f[0], f[1], int(f[2]), f[3], f[4], float(f[9])
        out[(chrom.removeprefix("chr"), pos, ref, alt)] = (gene, score)
    return out


def clinvar_labels(genes_tsv: Path) -> dict[tuple[str, int, str, str], str]:
    labels = {}
    for line in genes_tsv.read_text().splitlines():
        chrom, start, end, _gene = line.split("\t")
        out = subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/MC\n",
                str(CLINVAR),
                "-r",
                f"{chrom}:{start}-{end}",
            ],
            capture_output=True,
            text=True,
            check=False,
        ).stdout
        for row in out.splitlines():
            c, p, ref, alt, sig, mc = row.split("\t")
            if "missense" not in mc or "Conflicting" in sig:
                continue
            if "athogenic" in sig:
                labels[(c, int(p), ref, alt)] = "P"
            elif sig in ("Benign", "Likely_benign", "Benign/Likely_benign"):
                labels[(c, int(p), ref, alt)] = "B"
    return labels


def lr_at(threshold: float, pos: list[float], neg: list[float]) -> tuple[float, str, int, int]:
    """Positive likelihood ratio at ``threshold``, with a Rule-of-Three floor."""
    p_hits = sum(1 for s in pos if s >= threshold)
    n_hits = sum(1 for s in neg if s >= threshold)
    sens = p_hits / len(pos)
    if n_hits == 0:
        # 0 of n benign clear the bar: the rate is below 3/n with ~95% confidence.
        fpr = 3.0 / len(neg)
        bounded = True
    else:
        fpr = n_hits / len(neg)
        bounded = False
    lr = sens / fpr if fpr else float("inf")
    return lr, ("<=" if bounded else "="), p_hits, n_hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", type=Path, required=True, help="gene<TAB>AlphaMissense rows for the family")
    ap.add_argument("--genes", type=Path, required=True, help="chrom<TAB>start<TAB>end<TAB>gene")
    ap.add_argument("--query", default="15:40220612:T:G", help="the variant to place, chrom:pos:ref:alt")
    ap.add_argument("--bootstrap", type=int, default=2000)
    args = ap.parse_args()

    scores = load_scores(args.scores)
    labels = clinvar_labels(args.genes)

    pos = [scores[k][1] for k, v in labels.items() if v == "P" and k in scores]
    neg = [scores[k][1] for k, v in labels.items() if v == "B" and k in scores]
    print("=" * 78)
    print("Family calibration — spindle assembly checkpoint / kinetochore genes")
    print("=" * 78)
    print(f"  genes                    : {len(args.genes.read_text().splitlines())}")
    print(f"  AlphaMissense scores      : {len(scores):,}")
    print(f"  ClinVar P/LP missense     : {len(pos)} scored")
    print(f"  ClinVar B/LB missense     : {len(neg)} scored")
    if len(pos) < 5 or len(neg) < 20:
        print("\n  🔴 too few labelled variants to calibrate; reporting distributions only")

    print(f"\n  P/LP  median {statistics.median(pos):.4f}   min {min(pos):.4f}  max {max(pos):.4f}")
    print(f"  B/LB  median {statistics.median(neg):.4f}   min {min(neg):.4f}  max {max(neg):.4f}")

    # Separation: probability a random pathogenic outscores a random benign.
    wins = sum((a > b) + 0.5 * (a == b) for a in pos for b in neg)
    auc = wins / (len(pos) * len(neg))
    print(f"  separation (AUC)          : {auc:.3f}")

    chrom, p, ref, alt = args.query.split(":")
    key = (chrom, int(p), ref, alt)
    if key not in scores:
        print(f"\n  query {args.query} has no AlphaMissense score")
        return 1
    gene, qscore = scores[key]
    print(f"\n  query {args.query} ({gene}) AlphaMissense = {qscore:.4f}")

    print("\n" + "=" * 78)
    print("Likelihood ratio at the query's score, and at the generic boundary")
    print("=" * 78)
    print(f"{'threshold':>12} {'sens':>7} {'benign>=t':>10} {'LR+':>12}   evidence strength")
    for label, t in (("query score", qscore), ("AlphaMissense 0.564", 0.564), ("0.90", 0.90), ("0.80", 0.80)):
        lr, rel, ph, nh = lr_at(t, pos, neg)
        print(f"{t:>12.4f} {ph}/{len(pos):<5} {nh:>10} {rel}{lr:>10.1f}   {strength_for(lr):<18} [{label}]")

    # ---- bootstrap over the positives -------------------------------------
    rng = random.Random(0)
    lrs = []
    for _ in range(args.bootstrap):
        sample = [rng.choice(pos) for _ in pos]
        lrs.append(lr_at(qscore, sample, neg)[0])
    lrs.sort()
    lo, hi = lrs[int(0.025 * len(lrs))], lrs[int(0.975 * len(lrs))]
    print(f"\n  bootstrap over {len(pos)} positives ({args.bootstrap} resamples):")
    print(f"    LR+ 95% interval  {lo:.1f} to {hi:.1f}")
    print(f"    strength at the lower bound: {strength_for(lo)}")

    print("\n" + "=" * 78)
    print("🔴 What this does NOT license")
    print("=" * 78)
    print("  This calibration was built AFTER seeing the query's score. ClinGen requires")
    print("  the tool and threshold to be pre-specified, so it cannot be used to claim")
    print("  PP3 for this variant, and the submission does not claim it. What it does")
    print("  show is how far the generic genome-wide thresholds sit from what this gene")
    print("  family's own labelled data implies — which is the gap a VCEP would close.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
