#!/usr/bin/env python3
"""Screen for mosaic whole-chromosome aneuploidy from a single-sample VCF.

The defining phenotype of MVA is *variegated mosaic aneuploidy*, normally
measured by karyotyping cultured cells. This asks whether it is visible in the
sequencing data itself, from two signals available in a variants-only VCF:

* **B-allele frequency at heterozygous sites.** In a diploid region het BAF
  centres on 0.5. A whole-chromosome gain present in a fraction *f* of cells
  splits it into two clusters at ``(1+f)/(2+f)`` and ``1/(2+f)``.
* **Read depth at those sites**, relative to the genome median.

🔴 **Read the sign of the shift before believing it.** Short-read alignment has a
systematic preference for the reference allele, which pulls *every* chromosome's
mean BAF slightly below 0.5. A genuine mosaic gain is **bimodal** and affects one
or a few chromosomes; reference bias is **one-sided and genome-wide**. A run in
which all shifts are negative is measuring the aligner, not the karyotype.

Memory: per-chromosome running sums only. Nothing is held genome-wide.
"""

from __future__ import annotations

import argparse
import collections
import math
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

CANONICAL = [str(i) for i in range(1, 23)] + ["X", "Y"]
#: Quality floor for a site to count. Deliberately strict: this measurement is
#: a small shift in a mean, so noisy sites cost more than missing sites.
INCLUDE = 'TYPE="snp" && FILTER="PASS" && GT="het" && FMT/DP>=25 && FMT/GQ>=50'


def implied_mosaic_fraction(shift: float) -> float:
    """Cell fraction of a whole-chromosome gain that would produce ``shift``.

    A gain in fraction *f* of cells moves the gained allele's BAF from 0.5 to
    ``(1+f)/(2+f)``, i.e. a shift of ``f / (2(2+f))``. Inverted here.
    """
    d = abs(shift)
    return (2 * d) / (0.5 - d) if d < 0.5 else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", type=Path, default=REPO / "data/private/raw/proband.vcf.gz")
    ap.add_argument("--min-sites", type=int, default=1000)
    args = ap.parse_args()

    if not args.vcf.exists():
        print(f"VCF not found: {args.vcf}", file=sys.stderr)
        return 1

    # chrom -> [n, sum(baf), sum(baf^2), sum(dp)]
    stats: dict[str, list[float]] = collections.defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
    proc = subprocess.Popen(
        ["bcftools", "query", "-i", INCLUDE, "-f", "%CHROM\t[%AD]\t[%DP]\n", str(args.vcf)],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        try:
            chrom, ad, dp = line.rstrip("\n").split("\t")
        except ValueError:
            continue
        if chrom not in CANONICAL:
            continue
        parts = ad.split(",")
        try:
            ref, alt, depth = int(parts[0]), int(parts[1]), int(dp)
        except (ValueError, IndexError):
            continue
        total = ref + alt
        if total < 25:
            continue
        baf = alt / total
        s = stats[chrom]
        s[0] += 1
        s[1] += baf
        s[2] += baf * baf
        s[3] += depth
    proc.wait()

    total_sites = int(sum(s[0] for s in stats.values()))
    genome_depth = sum(s[3] for s in stats.values()) / max(total_sites, 1)
    print(f"{total_sites:,} high-confidence het SNVs; mean depth at het sites {genome_depth:.1f}\n")

    header = (
        f"{'chr':>4} {'n_het':>9} {'meanBAF':>9} {'SEM':>8} {'shift':>9} {'DPratio':>8} {'implied f':>10}"
    )
    print(header)
    print("-" * len(header))

    autosomal_means: list[float] = []
    for chrom in CANONICAL:
        if chrom not in stats:
            continue
        n, s1, s2, s3 = stats[chrom]
        if n < args.min_sites:
            continue
        mean = s1 / n
        var = max(s2 / n - mean * mean, 0.0)
        sem = math.sqrt(var / n)
        depth_ratio = (s3 / n) / genome_depth
        if chrom.isdigit():
            autosomal_means.append(mean)
        print(
            f"{chrom:>4} {int(n):>9,} {mean:>9.5f} {sem:>8.5f} {mean - 0.5:>+9.5f} "
            f"{depth_ratio:>8.3f} {implied_mosaic_fraction(mean - 0.5):>10.3f}"
        )

    if len(autosomal_means) < 2:
        return 0
    grand = statistics.mean(autosomal_means)
    scatter = statistics.pstdev(autosomal_means)
    limit = implied_mosaic_fraction(3 * scatter)
    n_negative = sum(1 for m in autosomal_means if m < 0.5)

    print(f"\nautosomal mean BAF {grand:.5f}, between-chromosome sd {scatter:.5f}")
    print(f"detection limit (3 sd): mosaic whole-chromosome gain in ~{limit * 100:.0f}% of cells")
    if n_negative == len(autosomal_means):
        print(
            f"\n🔴 all {n_negative} autosomal shifts are NEGATIVE — one-sided and genome-wide.\n"
            "   That is reference-allele mapping bias, not aneuploidy. A real mosaic gain is\n"
            "   bimodal and chromosome-specific. Report this as a null result."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
