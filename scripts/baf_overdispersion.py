"""Test for *variegated* aneuploidy via BAF overdispersion.

The repo's existing screen tests each chromosome for a *shift* in mean BAF.
That is the right test for a clonal mosaic trisomy, and the wrong one for MVA:
in variegated aneuploidy every chromosome is gained in some cells and lost in
others, so the means stay at 0.5 by symmetry and only the *variance* moves.

For a het site at depth D in a pure diploid, BAF ~ Binomial(D, p)/D, so
Var = p(1-p)/D. Aneuploid cells at fraction f pull individual sites off 0.5,
adding variance on top. The ratio of observed to binomial-expected variance is
therefore the statistic with power against MVA; the per-chromosome mean is not.
"""

from __future__ import annotations

import math
import subprocess
from collections import defaultdict

VCF = "data/private/raw/proband.vcf.gz"
AUTOSOMES = [str(c) for c in range(1, 23)]
MIN_DP, MAX_DP, MIN_MQ = 20, 80, 55

cmd = [
    "bcftools",
    "query",
    "-f",
    "%CHROM\t%MQ\t[%GT\t%AD\t%DP\t%GQ]\n",
    "-i",
    f'TYPE="snp" && N_ALT=1 && FILTER="PASS" && GT="het" && MQ>={MIN_MQ} && GQ=99 '
    f"&& FMT/DP>={MIN_DP} && FMT/DP<={MAX_DP}",
    VCF,
]

per_chrom: dict[str, list[tuple[float, int]]] = defaultdict(list)
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
assert proc.stdout is not None
for line in proc.stdout:
    chrom, _mq, _gt, ad, dp, _gq = line.rstrip("\n").split("\t")
    if chrom not in AUTOSOMES:
        continue
    try:
        ref_c, alt_c = (int(x) for x in ad.split(",")[:2])
    except ValueError:
        continue
    tot = ref_c + alt_c
    if tot < MIN_DP or tot != int(dp):
        continue
    per_chrom[chrom].append((alt_c / tot, tot))
proc.wait()

print(
    f"{'chrom':>6} {'n_het':>8} {'meanBAF':>8} {'varObs':>9} {'varBinom':>9} {'ratio':>7} {'implied f':>10}"
)
all_obs, all_exp, all_n = 0.0, 0.0, 0
rows = []
for c in AUTOSOMES:
    v = per_chrom.get(c, [])
    if len(v) < 1000:
        continue
    n = len(v)
    mean = sum(b for b, _ in v) / n
    var_obs = sum((b - mean) ** 2 for b, _ in v) / (n - 1)
    # Binomial expectation uses the observed mean, so reference-mapping bias
    # (which shifts the mean) does not leak into the overdispersion estimate.
    var_bin = sum(mean * (1 - mean) / d for _, d in v) / n
    excess = max(var_obs - var_bin, 0.0)
    # A whole-chromosome gain/loss in a fraction f of cells displaces BAF by
    # about f/4 in either direction, so excess SD ~= f/4 -> f ~= 4*sqrt(excess).
    implied_f = 4 * math.sqrt(excess)
    rows.append((c, n, mean, var_obs, var_bin, var_obs / var_bin, implied_f))
    all_obs += var_obs * n
    all_exp += var_bin * n
    all_n += n
    print(
        f"{c:>6} {n:>8} {mean:>8.4f} {var_obs:>9.5f} {var_bin:>9.5f} {var_obs / var_bin:>7.3f} {implied_f:>10.3f}"
    )

print(f"\ngenome-wide autosomal: n={all_n}  overdispersion ratio = {all_obs / all_exp:.3f}")
ratios = [r[5] for r in rows]
m = sum(ratios) / len(ratios)
sd = (sum((x - m) ** 2 for x in ratios) / (len(ratios) - 1)) ** 0.5
print(f"per-chromosome ratio: mean {m:.3f}  sd {sd:.3f}  min {min(ratios):.3f}  max {max(ratios):.3f}")
worst = max(rows, key=lambda r: r[5])
print(f"most overdispersed chromosome: {worst[0]}  ratio {worst[5]:.3f}  implied f {worst[6]:.3f}")
