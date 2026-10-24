"""Is the BUB1B gene body inside a run of homozygosity?

It matters: an ROH would mean consanguinity, which raises the prior on a
HOMOZYGOUS causal variant and would undercut a compound-heterozygous call.
The gene body carries far fewer variants than the chr15 average, so this has
to be excluded rather than assumed.
"""

from __future__ import annotations

import subprocess
from collections import defaultdict

VCF = "data/private/raw/proband.vcf.gz"
WIN = 60_000  # one BUB1B gene span
BUB1B = (40160984, 40221137)

out = subprocess.run(
    ["bcftools", "query", "-r", "15", "-f", "%POS\t[%GT]\n", "-i", 'FILTER="PASS"', VCF],
    capture_output=True,
    text=True,
    check=True,
).stdout

tot = defaultdict(int)
het = defaultdict(int)
for line in out.splitlines():
    pos, gt = line.split("\t")
    w = int(pos) // WIN
    tot[w] += 1
    if gt in ("0/1", "0|1", "1|0"):
        het[w] += 1

# Only windows with any callable content; chr15 is acrocentric so the p-arm is empty.
wins = [w for w in sorted(tot) if tot[w] > 0]
counts = sorted(tot[w] for w in wins)
hets = sorted(het[w] for w in wins)
n = len(counts)


def pct(sorted_vals, v):
    return 100.0 * sum(1 for x in sorted_vals if x < v) / len(sorted_vals)


bw = BUB1B[0] // WIN
print(f"chr15, {WIN // 1000} kb windows with calls: n={n}")
print(f"  median variants/window : {counts[n // 2]}")
print(f"  median hets/window     : {hets[n // 2]}")
print(f"  10th pct variants      : {counts[n // 10]}   90th pct: {counts[9 * n // 10]}")

print(f"\nBUB1B window ({bw * WIN}-{(bw + 1) * WIN}):")
print(f"  variants {tot[bw]:>4}  -> {pct(counts, tot[bw]):.1f}th percentile of chr15")
print(f"  hets     {het[bw]:>4}  -> {pct(hets, het[bw]):.1f}th percentile of chr15")

print("\nneighbouring windows (an ROH is a long run of het-poor windows):")
for w in range(bw - 6, bw + 7):
    if w not in tot:
        continue
    mark = "  <== BUB1B" if w == bw else ""
    bar = "#" * min(40, het[w] // 2)
    print(f"  {w * WIN / 1e6:7.2f} Mb  var {tot[w]:>4}  het {het[w]:>4} {bar}{mark}")

# ROH definition: a run of windows with het density < 10% of the chromosome median.
thresh = max(1, hets[n // 2] // 10)
runs, cur = [], 0
for w in wins:
    if het[w] < thresh:
        cur += 1
    else:
        if cur:
            runs.append(cur)
        cur = 0
if cur:
    runs.append(cur)
print(f"\nROH threshold: <{thresh} hets per {WIN // 1000} kb window")
print(
    f"longest het-poor run on chr15: {max(runs) if runs else 0} windows "
    f"({(max(runs) if runs else 0) * WIN / 1e6:.2f} Mb)"
)
print(
    f"BUB1B window het count {het[bw]} vs threshold {thresh} -> "
    f"{'INSIDE an ROH' if het[bw] < thresh else 'NOT in an ROH'}"
)
