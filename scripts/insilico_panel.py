#!/usr/bin/env python3
"""Score BUB1B missense alleles with the ClinGen-calibrated predictors.

Why this exists
---------------
A missense predictor's score is only useful as ACMG evidence if the score has
been calibrated to an evidence strength. ClinGen's Sequence Variant
Interpretation group did that for 13 tools (Pejaver et al., *Am J Hum Genet*
109:2163-2177, 2022, Table 2) and later for AlphaMissense, ESM1b and VARITY
(Bergquist et al., *Genet Med* 2025, Table 1). This script scores the child's
missense allele with REVEL, VEST4, PolyPhen-2 and SIFT, and reports the strength
each tool's score would license. AlphaMissense is handled in
``alphamissense_bub1b.py``.

The script reports every tool. It does not pick the one that agrees with us:
ClinGen asks for a single tool to be chosen before the variant is looked at.

Everything is read from **bulk public files downloaded in full and queried
locally**. No proband-derived coordinate is transmitted. Fetch them with
``make insilico-refs``.

Controls matter more than the headline
--------------------------------------
A predictor agreeing with us proves nothing on its own. The script scores:

* ``p.Leu1012Pro`` — ClinVar variation 6765, the canonical *BUB1B* kinase-domain
  missense of the MVA1 literature, **10 residues** from our allele. A positive
  control drawn from the disease, not from a benchmark set.
* ``p.Gln467His`` — the *only* ClinVar Likely-pathogenic missense in the whole
  gene. A second positive control, and the one our tools disagree about.
* every ClinVar Benign/Likely-benign missense in the gene — negative controls,
  which is where a predictor that simply scores everything highly gets caught.
"""

from __future__ import annotations

import argparse
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REFS = REPO / "data/reference/insilico"
CLINVAR = REPO / "data/reference/clinvar/clinvar.vcf.gz"

#: BUB1B, GRCh38, from the Ensembl GFF3 — a public gene boundary, not a
#: proband-derived coordinate.
CONTIG, GENE_START, GENE_END = "15", 40160984, 40221137
MANE = "ENST00000287598"


#: Pejaver et al. 2022, Table 2 ("Estimated threshold ranges for all tools in
#: this study corresponding to the four pathogenic and four benign intervals",
#: p. 2169). Transcribed from the paper, not from memory. ``None`` where the
#: tool did not reach that strength.
#:
#: Each entry: (PP3_supporting, PP3_moderate, PP3_strong) as lower bounds, and
#: (BP4_supporting, BP4_moderate, BP4_strong, BP4_verystrong) as upper bounds.
#: ``higher_is_worse`` is False for SIFT, where a low score means damaging.
@dataclass(frozen=True)
class Calibration:
    name: str
    pp3: tuple[float | None, float | None, float | None]
    bp4: tuple[float | None, float | None, float | None, float | None]
    higher_is_worse: bool = True

    def strength(self, score: float | None) -> str:
        """Return the ACMG criterion and strength this score licenses."""
        if score is None:
            return "no score"
        cmp = (lambda a, b: a >= b) if self.higher_is_worse else (lambda a, b: a <= b)
        for bound, label in zip(
            reversed(self.pp3), ("PP3_strong", "PP3_moderate", "PP3_supporting"), strict=True
        ):
            if bound is not None and cmp(score, bound):
                return label
        rev = (lambda a, b: a <= b) if self.higher_is_worse else (lambda a, b: a >= b)
        for bound, label in zip(
            self.bp4, ("BP4_supporting", "BP4_moderate", "BP4_strong", "BP4_verystrong"), strict=True
        ):
            if bound is not None and rev(score, bound):
                return label
        return "indeterminate"


CALIBRATED = {
    "REVEL": Calibration("REVEL", (0.644, 0.773, 0.932), (0.290, 0.183, 0.016, 0.003)),
    "VEST4": Calibration("VEST4", (0.764, 0.861, 0.965), (0.449, 0.302, None, None)),
    "PolyPhen2": Calibration("PolyPhen2", (0.978, 0.999, None), (0.113, 0.009, None, None)),
    "SIFT": Calibration("SIFT", (0.001, 0.0, None), (0.327, None, None, None), higher_is_worse=False),
}

#: AlphaMissense's own class boundaries are 0.34 and 0.564. They are not ACMG
#: evidence strengths. The ClinGen-calibrated bands (Bergquist 2025, Table 1) start
#: at 0.792 for PP3_supporting; see ``alphamissense_bub1b.py``.
AM_LIKELY_PATH = 0.564


# ---------------------------------------------------------------------------
# Loaders — one per bulk file format
# ---------------------------------------------------------------------------


def load_alphamissense(path: Path) -> dict[tuple[int, str, str], tuple[str, float, str]]:
    out = {}
    for line in path.read_text().splitlines():
        _c, pos, ref, alt, _g, _u, _t, pv, score, cls = line.split("\t")
        out[(int(pos), ref, alt)] = (pv, float(score), cls)
    return out


def load_revel(path: Path) -> dict[tuple[int, str, str], float]:
    out: dict[tuple[int, str, str], float] = {}
    for line in path.read_text().splitlines():
        f = line.split(",")
        if f[2] == "." or not f[7]:
            continue
        key = (int(f[2]), f[3], f[4])
        out[key] = max(out.get(key, 0.0), float(f[7]))
    return out


def _bed_groups(fields: list[str], start: int, width: int) -> list[list[str]]:
    rest = fields[start:]
    return [rest[i : i + width] for i in range(0, len(rest) - width + 1, width)]


def load_vest(path: Path) -> dict[tuple[int, str, str], float]:
    """UCSC dbNsfpVest.bb: ref, ensembl-ids, all-ids, then (alt, aa, scores)*."""
    out: dict[tuple[int, str, str], float] = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        pos, ref, ids = int(f[2]), f[3], f[5].split(",")
        idx = ids.index(MANE) if MANE in ids else None
        for alt, _aa, scores in _bed_groups(f, 6, 3):
            vals = [s for s in scores.split(",") if s not in ("", ".")]
            if not vals:
                continue
            # Prefer the MANE transcript; fall back to the most damaging call.
            pick = vals[idx] if idx is not None and idx < len(vals) else max(vals, key=float)
            out[(pos, ref, alt)] = float(pick)
    return out


def load_polyphen2(path: Path) -> dict[tuple[int, str, str], float]:
    """UCSC dbNsfpPolyPhen2.bb: (alt, HDIV, HDIVpred, HVAR, HVARpred)*.

    Pejaver et al. calibrated PolyPhen-2 **HumVar**, so HVAR is what is read.
    """
    out: dict[tuple[int, str, str], float] = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        pos, ref = int(f[2]), f[3]
        for alt, _hdiv, _hdp, hvar, _hvp in _bed_groups(f, 5, 5):
            # PolyPhen-2 reports one score per protein isoform, comma-joined.
            # There is no transcript key on this track, so take the most
            # damaging - the conservative direction for a negative control.
            vals = [v for v in hvar.split(",") if v not in ("", ".")]
            if not vals:
                continue
            out[(pos, ref, alt)] = max(out.get((pos, ref, alt), 0.0), max(map(float, vals)))
    return out


def load_sift(path: Path) -> dict[tuple[int, str, str], float]:
    """UCSC dbNsfpSift.bb: ref, ensembl-ids, then (alt, scores)*. Low = damaging."""
    out: dict[tuple[int, str, str], float] = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        pos, ref, ids = int(f[2]), f[3], f[4].split(",")
        idx = ids.index(MANE) if MANE in ids else None
        for alt, scores in _bed_groups(f, 5, 2):
            vals = [s for s in scores.split(",") if s not in ("", ".")]
            if not vals:
                continue
            pick = vals[idx] if idx is not None and idx < len(vals) else min(vals, key=float)
            out[(pos, ref, alt)] = float(pick)
    return out


def clinvar_missense() -> list[tuple[int, str, str, str, str]]:
    """Every ClinVar missense record in the gene: (pos, ref, alt, id, clnsig)."""
    out = subprocess.run(
        [
            "bcftools",
            "query",
            "-f",
            "%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/MC\n",
            str(CLINVAR),
            "-r",
            f"{CONTIG}:{GENE_START}-{GENE_END}",
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    rows = []
    for line in out.splitlines():
        pos, ref, alt, vid, sig, mc = line.split("\t")
        if "missense" in mc:
            rows.append((int(pos), ref, alt, vid, sig))
    return rows


# ---------------------------------------------------------------------------

WATCHLIST = [
    ((40220612, "T", "G"), "p.Asn1002Lys", "OUR ALLELE 2 — absent from ClinVar as T>G"),
    ((40220612, "T", "A"), "p.Asn1002Lys", "same protein change, ClinVar 4600147 = VUS"),
    ((40220641, "T", "C"), "p.Leu1012Pro", "ClinVar 6765 — canonical MVA1 kinase allele"),
    ((40199727, "G", "T"), "p.Gln467His", "the only ClinVar Likely-pathogenic missense"),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=REFS, help="directory of BUB1B predictor slices")
    args = ap.parse_args()

    am = load_alphamissense(args.refs / "bub1b_alphamissense.tsv")
    scores = {
        "REVEL": load_revel(args.refs / "bub1b_revel.csv"),
        "VEST4": load_vest(args.refs / "bub1b_Vest.bed"),
        "PolyPhen2": load_polyphen2(args.refs / "bub1b_PolyPhen2.bed"),
        "SIFT": load_sift(args.refs / "bub1b_Sift.bed"),
    }

    print("=" * 100)
    print("BUB1B missense — ClinGen SVI calibrated predictors (Pejaver 2022 thresholds)")
    print("=" * 100)
    print(f"{'variant':<14} {'protein':<13} {'REVEL':>21} {'VEST4':>21} {'PolyPhen2':>19} {'SIFT':>19}")
    for key, prot, _note in WATCHLIST:
        row = f"{CONTIG}:{key[0]} {key[1]}>{key[2]}"
        cells = []
        for tool, cal in CALIBRATED.items():
            s = scores[tool].get(key)
            cells.append(f"{'--' if s is None else f'{s:.3f}'} {cal.strength(s):<15}")
        print(
            f"{row:<14} {prot:<13} "
            + " ".join(f"{c:>21}" for c in cells[:2])
            + " "
            + " ".join(f"{c:>19}" for c in cells[2:])
        )

    print()
    print("Where each allele sits in that tool's own BUB1B distribution (percentile,")
    print("100 = most damaging by that tool). A tool can rank a variant highly and")
    print("still decline PP3: the ClinGen threshold is absolute, not relative.")
    print(f"{'variant':<14} {'protein':<13} " + " ".join(f"{k:>11}" for k in CALIBRATED))
    for key, prot, _note in WATCHLIST:
        cells = []
        for tool, cal in CALIBRATED.items():
            s = scores[tool].get(key)
            if s is None:
                cells.append("--")
                continue
            vals = sorted(scores[tool].values())
            worse = sum(1 for v in vals if (v < s if cal.higher_is_worse else v > s))
            cells.append(f"{100.0 * worse / len(vals):.1f}")
        print(f"{CONTIG}:{key[0]} {key[1]}>{key[2]:<3} {prot:<13} " + " ".join(f"{c:>11}" for c in cells))

    print()
    print(
        "AlphaMissense (its own class labels below are not ACMG strengths; the\n  ClinGen-calibrated bands are applied in alphamissense_bub1b.py):"
    )
    ranked = sorted(v[1] for v in am.values())
    for key, prot, note in WATCHLIST:
        hit = am.get(key)
        if not hit:
            continue
        pct = 100.0 * sum(1 for s in ranked if s < hit[1]) / len(ranked)
        print(
            f"  {CONTIG}:{key[0]} {key[1]}>{key[2]}  {prot:<13} {hit[1]:.4f}  {hit[2]:<18} {pct:5.1f}th pct   {note}"
        )

    # ---- controls ---------------------------------------------------------
    cv = clinvar_missense()
    benign = [(p, r, a) for p, r, a, _i, s in cv if "enign" in s and "Conflicting" not in s]
    print()
    print("=" * 100)
    print(f"NEGATIVE CONTROLS — every ClinVar Benign/Likely-benign missense in BUB1B (n={len(benign)})")
    print("=" * 100)
    print(f"{'tool':<12} {'n scored':>9} {'median':>9} {'max':>9}   {'called PP3 at any strength':>28}")
    for tool, cal in CALIBRATED.items():
        vals = [scores[tool][k] for k in benign if k in scores[tool]]
        if not vals:
            continue
        bad = [v for v in vals if cal.strength(v).startswith("PP3")]
        print(
            f"{tool:<12} {len(vals):>9} {statistics.median(vals):>9.3f} "
            f"{(max(vals) if cal.higher_is_worse else min(vals)):>9.3f}   {len(bad):>10} of {len(vals):<15}"
        )
    am_vals = [am[k][1] for k in benign if k in am]
    am_bad = [v for v in am_vals if v > AM_LIKELY_PATH]
    print(
        f"{'AlphaMissense':<12} {len(am_vals):>9} {statistics.median(am_vals):>9.3f} {max(am_vals):>9.3f}   "
        f"{len(am_bad):>10} of {len(am_vals):<15} (>0.564, its own class boundary)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
