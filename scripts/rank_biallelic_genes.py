#!/usr/bin/env python3
"""Rank genes by biallelic-compatible rare damaging genotype — phenotype-blind.

Consumes the TSV written by ``scripts/coding_consequence_scan.sh`` and asks one
question, with no knowledge of the proband's phenotype and no candidate gene
supplied: *which genes carry a recessive genotype that could explain a severe
paediatric presentation?*

🔴 Three filters do almost all of the work, and leaving any of them out makes
the ranking useless rather than merely noisier. Run without them, the top of the
list is ``FCGBP``, ``MUC3A``, ``MUC4``, ``HLA-*``, ``LILRB*`` and the ``KRTAP``
cluster — polymorphic, segmentally duplicated, mismapping-prone genes in which
apparent homozygous loss-of-function is an alignment artefact, not biology:

1. **Mapping quality.** ``MQ >= 55``. Multi-mapping regions manufacture
   false homozygous LoF.
2. **Depth sanity.** ``15 <= DP <= 100``. Collapsed repeats show inflated depth.
3. **Rarity.** A real gnomAD allele frequency where we hold one, and the
   "no dbSNP rsID in the proband VCF" proxy elsewhere.

🔴 **The rsID proxy is much weaker than it looks, and it was measured, not
assumed.** Audited against gnomAD v4.1 exomes on chr15, only **35%** of the
variants it calls rare actually are (6 of 17; 11 are commoner than 0.1%) — and
that single-chromosome figure is optimistic: measured genome-wide with
``--gnomad-tsv`` the proxy agrees only **26.8%** of the time (227 of 847).
Prefer ``--gnomad-tsv``; the proxy is a fallback for when no frequency table is
available, not a frequency filter. It
fails in the safe direction — it over-retains rather than over-excludes, so it
never hides the answer — but it is a screen, not a frequency filter, and the
ranking it produces is correspondingly noisy. Where gnomAD is available the real
AF replaces it, and chr15 then collapses from 480 damaging variants to 9, and
from 15 candidate genes to 2.

🔴 **"Absent from gnomAD" does not mean "rare".** It usually means the allele is
written differently on the two sides — an unnormalised indel, or a multi-allelic
site. Defaulting absence to *keep* once let an unnormalised indel with a gnomAD AF
of **0.52** through as a "rare homozygous LoF" and rank it above the real
answer. Absence is only treated as rarity when the proband VCF also gives
the site no rsID; a known rsID with no gnomAD match is a representation
mismatch, not a discovery. Normalise with ``bcftools norm`` before trusting any
cross-database comparison.

Pathogenicity evidence comes from a **local** ClinVar VCF and frequency from a
**local** gnomAD VCF. Nothing in this script touches the network, so no
proband-derived coordinate leaves the machine.
"""

from __future__ import annotations

import argparse
import collections
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

LOF = frozenset({"stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost", "stop_lost"})
#: Most-severe-first. Index position is the severity rank.
SEVERITY: tuple[str, ...] = (
    "stop_gained",
    "frameshift",
    "splice_acceptor",
    "splice_donor",
    "start_lost",
    "stop_lost",
    "inframe_deletion",
    "inframe_insertion",
    "missense",
)

#: Biallelic genotype patterns, most convincing first.
PATTERNS: tuple[tuple[str, int, str], ...] = (
    ("hom-LoF", 100, "homozygous loss of function"),
    ("2x het-LoF", 95, "two heterozygous LoF alleles"),
    ("het-LoF + het-missense", 85, "one LoF allele plus one missense allele"),
    ("hom-missense", 60, "homozygous missense"),
    ("2x het-missense", 40, "two heterozygous missense alleles"),
)
CLINVAR_BONUS = 30
#: More than this many rare damaging alleles in one gene reads as artefact.
CROWDING_LIMIT = 4
CROWDING_PENALTY = 0.5


def most_severe(bcsq: str) -> tuple[str | None, str | None]:
    """Return the most severe consequence across transcripts, and its gene.

    An annotation carried only by an NMD-biotype transcript is not evidence that
    the protein changes: every isoform bearing it is non-productive. Such a
    variant is not protein-altering and is dropped.

    🔴 This used to be written as ``if "NMD_transcript" in consequence and best
    is not None: continue``, which made the outcome depend on the ORDER
    ``bcftools csq`` happened to emit transcripts in: an NMD annotation arriving
    first was accepted, because ``best`` was still None. That is not
    hypothetical — it put ``TP53`` into the genome-wide damaging set as a
    "frameshift" on the strength of ``ENST00000635293``, an NMD-biotype
    transcript, when no protein-coding transcript is affected at all.
    """
    best: str | None = None
    gene: str | None = None
    for entry in bcsq.split(","):
        parts = entry.split("|")
        if len(parts) < 2:
            continue
        consequence, symbol = parts[0], parts[1]
        if "NMD_transcript" in consequence:
            # Carried only by a non-productive isoform: the protein is unchanged.
            continue
        for candidate in SEVERITY:
            if candidate in consequence:
                if best is None or SEVERITY.index(candidate) < SEVERITY.index(best):
                    best, gene = candidate, symbol
                break
    return best, gene


def zygosity(genotype: str) -> str:
    alleles = genotype.replace("|", "/").split("/")
    if alleles == ["1", "1"]:
        return "hom"
    if "1" in alleles and "0" in alleles:
        return "het"
    return "other"


def load_variants(tsv: Path) -> list[dict]:
    out: list[dict] = []
    for line in tsv.open():
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        chrom, pos, ref, alt, _qual, gt, dp, ad, gq, bcsq = f[:10]
        consequence, gene = most_severe(bcsq)
        if not gene or not consequence:
            continue
        try:
            depth, quality = int(dp), int(gq)
        except ValueError:
            continue
        out.append(
            {
                "key": (chrom, int(pos), ref, alt),
                "gene": gene,
                "csq": consequence,
                "zyg": zygosity(gt),
                "dp": depth,
                "gq": quality,
                "ad": ad,
            }
        )
    return out


def bcftools_lookup(vcf: Path, keys: list[tuple], fmt: str) -> dict[tuple, list[str]]:
    """Query `vcf` at exactly the given positions. One subprocess, streamed."""
    regions = REPO / "data" / "private" / "work" / ".regions.tmp"
    regions.parent.mkdir(parents=True, exist_ok=True)
    with regions.open("w") as handle:
        for chrom, pos, _ref, _alt in sorted(set(keys), key=lambda k: (k[0], k[1])):
            handle.write(f"{chrom}\t{pos}\n")
    proc = subprocess.run(
        ["bcftools", "query", "-R", str(regions), "-f", fmt, str(vcf)],
        capture_output=True,
        text=True,
        check=False,
    )
    result: dict[tuple, list[str]] = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 4:
            continue
        result[(f[0], int(f[1]), f[2], f[3])] = f[4:]
    regions.unlink(missing_ok=True)
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", type=Path, default=REPO / "data/private/work/damaging_genomewide.tsv")
    ap.add_argument("--vcf", type=Path, default=REPO / "data/private/raw/proband.vcf.gz")
    ap.add_argument("--clinvar", type=Path, default=REPO / "data/reference/clinvar/clinvar.vcf.gz")
    ap.add_argument(
        "--no-clinvar",
        action="store_true",
        help="rank on the genotype pattern alone, with every ClinVar Pathogenic/Likely-pathogenic "
        "flag withheld. This is the ablation behind the report's statement that no single "
        "evidence axis solves this case. Without prior disease evidence the answer drops from "
        "rank 1 to rank 41, which is worth being able to reproduce rather than assert.",
    )
    ap.add_argument(
        "--gnomad",
        type=Path,
        default=REPO / "data/reference/gnomad/gnomad.exomes.v4.1.sites.chr15.vcf.bgz",
        help="single-contig gnomAD sites VCF (chr15). Audits the rsID rarity proxy on one "
        "chromosome rather than replacing it. Superseded by --gnomad-tsv, which replaces "
        "the proxy genome-wide; this remains the low-setup path.",
    )
    ap.add_argument(
        "--gnomad-tsv",
        type=Path,
        default=None,
        help="precomputed genome-wide gnomAD AF table (CHROM POS REF ALT AF nhomalt AC AN), as "
        "written by scripts/gnomad_af_genomewide.sh (make gnomad-af). When given, a REAL allele "
        "frequency is used on every contig and the rsID rarity proxy is not used at all.",
    )
    ap.add_argument("--max-af", type=float, default=1e-3, help="rarity threshold for the audit")
    ap.add_argument("--min-mq", type=float, default=55.0)
    ap.add_argument(
        "--min-ab",
        type=float,
        default=0.25,
        help="minimum alternate-allele fraction for a HETEROZYGOUS call. A true germline "
        "het sits near 0.5; below ~0.25 a call is far more likely an alignment artefact. "
        "MQ and DP do not catch this: a call can pass DP 34 / GQ 99 / MQ 58 at an "
        "alternate fraction of 0.147 and still be an alignment artefact.",
    )
    ap.add_argument("--min-dp", type=int, default=15)
    ap.add_argument("--max-dp", type=int, default=100)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument(
        "--contig",
        default=None,
        help="restrict to one contig, e.g. 15. Use with a matching --gnomad file to see the "
        "ranking under a real allele-frequency filter rather than the rsID proxy.",
    )
    ap.add_argument("--highlight", default="BUB1B", help="report this gene's rank explicitly")
    args = ap.parse_args()

    variants = load_variants(args.tsv)
    if args.contig:
        variants = [v for v in variants if v["key"][0] == args.contig]
    print(
        f"loaded {len(variants):,} damaging variants" + (f" on contig {args.contig}" if args.contig else ""),
        file=sys.stderr,
    )

    keys = [v["key"] for v in variants]
    site = bcftools_lookup(args.vcf, keys, "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/MQ\n")
    if args.no_clinvar:
        # Ablation: withhold prior disease evidence entirely, so the ranking rests
        # on genotype pattern and rarity alone.
        clinvar = {}
        print("ClinVar evidence WITHHELD (--no-clinvar): ranking on genotype alone", file=sys.stderr)
    else:
        clinvar = bcftools_lookup(args.clinvar, keys, "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\n")

    # gnomAD v4 uses chr-prefixed contigs; the proband VCF does not.
    # Values are (max AF across ALT alleles, nhomalt). AF can be a comma list at a
    # multi-allelic site, so parse every field and keep the largest - silently
    # dropping an unparseable AF is what let a 52%-frequency allele through.
    gnomad: dict[tuple, tuple[float, str]] = {}
    gnomad_contigs: set[str] = set()
    if args.gnomad_tsv and args.gnomad_tsv.exists():
        # Genome-wide real frequencies: the rsID proxy is retired entirely.
        for line in args.gnomad_tsv.open():
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            contig = f[0].removeprefix("chr")
            values = []
            for piece in f[4].split(","):
                try:
                    values.append(float(piece))
                except ValueError:
                    continue
            if not values:
                continue
            gnomad_contigs.add(contig)
            gnomad[(contig, int(f[1]), f[2], f[3])] = (max(values), f[5])
        if not gnomad:
            raise SystemExit(
                f"{args.gnomad_tsv} parsed to 0 usable records. Refusing to continue: "
                "silently falling back to the rsID proxy here would look like a "
                "genome-wide real-frequency run and would not be one."
            )
        # 🔴 Retire the proxy on EVERY contig we are ranking, not merely on the
        # contigs that happen to have records. A contig with no gnomAD match is
        # then read as absence (rare, if the site also has no rsID), which is the
        # documented behaviour - rather than silently reverting to the proxy.
        gnomad_contigs = {k[0] for k in keys}
        print(
            f"loaded {len(gnomad):,} gnomAD records from {args.gnomad_tsv}; "
            f"real frequencies in force on all {len(gnomad_contigs)} contigs - rsID proxy DISABLED",
            file=sys.stderr,
        )
    elif args.gnomad.exists():
        gnomad_contigs = {"15"}
        wanted = [k for k in keys if k[0] in gnomad_contigs]
        raw = bcftools_lookup(
            args.gnomad,
            [(f"chr{c}", p, r, a) for c, p, r, a in wanted],
            "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/nhomalt\n",
        )
        for (contig, pos, ref, alt), fields in raw.items():
            values = []
            for piece in (fields[0] if fields else "").split(","):
                try:
                    values.append(float(piece))
                except ValueError:
                    continue
            if values:
                nhomalt = fields[1] if len(fields) > 1 else "?"
                gnomad[(contig.removeprefix("chr"), pos, ref, alt)] = (max(values), nhomalt)

    pathogenic = {k for k, v in clinvar.items() if v and "athogenic" in v[0] and "Conflicting" not in v[0]}

    kept: list[dict] = []
    proxy_rare: list[tuple] = []
    dropped = collections.Counter()
    for v in variants:
        meta = site.get(v["key"])
        if not meta:
            dropped["no-site-record"] += 1
            continue
        rsid = meta[0]
        try:
            mq = float(meta[1])
        except (IndexError, ValueError):
            mq = 0.0
        if mq < args.min_mq:
            dropped["low-MQ"] += 1
            continue
        if not (args.min_dp <= v["dp"] <= args.max_dp):
            dropped["depth-out-of-range"] += 1
            continue
        if v["zyg"] == "het":
            try:
                ref_c, alt_c = (int(x) for x in v["ad"].split(",")[:2])
            except ValueError:
                ref_c = alt_c = 0
            if ref_c + alt_c and alt_c / (ref_c + alt_c) < args.min_ab:
                dropped[f"allele-balance < {args.min_ab:g} (likely artefact)"] += 1
                continue

        v["af"] = v["nhomalt"] = None
        # Record what the rsID proxy *alone* would have decided, so the audit
        # below can score it against real frequencies rather than against the
        # post-AF-filter survivors (which would be circular).
        if v["key"][0] in gnomad_contigs and rsid == ".":
            proxy_rare.append(v["key"])
        if v["key"][0] in gnomad_contigs:
            # A real frequency is available for this contig, so use it.
            if v["key"] in gnomad:
                af, nhomalt = gnomad[v["key"]]
                if af > args.max_af:
                    dropped[f"common (gnomAD AF > {args.max_af:g})"] += 1
                    continue
                v["af"], v["nhomalt"] = af, nhomalt
            elif rsid != ".":
                # Known allele with no gnomAD match at this exact spelling: a
                # representation mismatch, not a rare variant. Do not keep it.
                dropped["no gnomAD match but has rsID (unnormalised)"] += 1
                continue
        elif rsid != ".":
            dropped["common (has rsID proxy)"] += 1
            continue

        v["clinvar_plp"] = v["key"] in pathogenic
        kept.append(v)

    print(
        f"kept {len(kept):,} rare, well-mapped damaging variants; dropped "
        + ", ".join(f"{n:,} {why}" for why, n in dropped.most_common()),
        file=sys.stderr,
    )

    # --- audit the rsID rarity proxy against real frequencies ---------------
    # Scores the proxy on its own terms: of the variants it alone would have
    # called rare (on a contig where we hold gnomAD), how many really are?
    if gnomad_contigs:
        seen = [k for k in proxy_rare if k in gnomad]
        rare = sum(1 for k in seen if gnomad[k][0] <= args.max_af)
        absent = len(proxy_rare) - len(seen)
        agree = absent + rare
        total = max(len(proxy_rare), 1)
        print(
            f"\nrsID-proxy audit on {'/'.join(sorted(gnomad_contigs))} "
            f"({len(proxy_rare)} variants the proxy alone calls rare):\n"
            f"  present, AF <= {args.max_af:g}        : {rare}\n"
            f"  present, COMMONER than {args.max_af:g} : {len(seen) - rare}\n"
            f"  absent from gnomAD entirely  : {absent}\n"
            f"  proxy agrees with real AF    : {agree}/{total}"
            f" ({100 * agree / total:.1f}%)",
            file=sys.stderr,
        )

    by_gene: dict[str, list[dict]] = collections.defaultdict(list)
    for v in kept:
        by_gene[v["gene"]].append(v)

    scored = []
    for gene, vs in by_gene.items():
        hom_lof = [v for v in vs if v["zyg"] == "hom" and v["csq"] in LOF]
        het_lof = [v for v in vs if v["zyg"] == "het" and v["csq"] in LOF]
        hom_mis = [v for v in vs if v["zyg"] == "hom" and v["csq"] not in LOF]
        het_mis = [v for v in vs if v["zyg"] == "het" and v["csq"] not in LOF]

        if hom_lof:
            pattern, base = PATTERNS[0][0], PATTERNS[0][1]
        elif len(het_lof) >= 2:
            pattern, base = PATTERNS[1][0], PATTERNS[1][1]
        elif len(het_lof) == 1 and het_mis:
            pattern, base = PATTERNS[2][0], PATTERNS[2][1]
        elif hom_mis:
            pattern, base = PATTERNS[3][0], PATTERNS[3][1]
        elif len(het_mis) >= 2:
            pattern, base = PATTERNS[4][0], PATTERNS[4][1]
        else:
            continue

        score = float(base)
        n_plp = sum(1 for v in vs if v["clinvar_plp"])
        if n_plp:
            score += CLINVAR_BONUS
        if len(vs) > CROWDING_LIMIT:
            score *= CROWDING_PENALTY
        scored.append((score, gene, pattern, n_plp, vs))

    scored.sort(key=lambda r: (-r[0], r[1]))
    print(f"\n{len(scored):,} genes carry a rare biallelic-compatible damaging genotype\n")
    print(f"{'rank':>5}  {'gene':<14} {'pattern':<24} {'alleles':>7} {'ClinVar P/LP':>13} {'score':>7}")
    print("-" * 78)
    for i, (score, gene, pattern, n_plp, vs) in enumerate(scored[: args.top], 1):
        print(f"{i:>5}  {gene:<14} {pattern:<24} {len(vs):>7} {n_plp:>13} {score:>7.1f}")

    for i, (score, gene, pattern, _n_plp, vs) in enumerate(scored, 1):
        if gene == args.highlight:
            print(f"\n>>> {gene}: rank {i} of {len(scored)}  [{pattern}, score {score:.1f}]")
            for v in sorted(vs, key=lambda v: v["key"][1]):
                c, p, r, a = v["key"]
                flags = " ClinVar-Pathogenic/LP" if v["clinvar_plp"] else ""
                if v["af"] is not None:
                    flags += f" gnomAD_AF={v['af']:.3g} nhomalt={v['nhomalt']}"
                elif c in gnomad_contigs:
                    flags += " absent-from-gnomAD (and no rsID)"
                print(f"      {c}:{p} {r}>{a}  {v['csq']:<12} {v['zyg']}  DP={v['dp']} AD={v['ad']}{flags}")
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
