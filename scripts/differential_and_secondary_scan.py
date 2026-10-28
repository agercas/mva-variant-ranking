#!/usr/bin/env python3
"""Exclude the competing diagnoses, and surface secondary findings — gene-agnostically.

Two questions, both answered without a gene panel:

**1. Does anything else in this genome explain the phenotype?**
Intersect every variant the proband carries against the *whole* of ClinVar's
Pathogenic / Likely-pathogenic set (~346,000 records), then ask which genes carry
a homozygous P/LP allele or two P/LP alleles — i.e. a complete recessive
genotype. If a competing diagnosis exists, it shows up here without anyone having
to think of it first.

**2. Are the phenotype-driven differentials actually excluded?**
For each named differential gene, print every damaging variant with its ClinVar
classification and rarity, so "excluded" is a statement with evidence under it
rather than an assertion.

🔴 Everything is local: ClinVar is a downloaded VCF, the proband VCF is in the
quarantine. No coordinate leaves the machine.
"""

from __future__ import annotations

import argparse
import collections
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

#: The genes linked to MVA, plus the reported genes. These are named in the results. The
#: other competing diagnoses are not listed here: they are read at run time from the
#: git-ignored EXTRA_DIFFERENTIAL_FILE (see load_differential), and the public view counts
#: them without naming them.
DIFFERENTIAL: dict[str, str] = {
    "LZTR1": "Noonan syndrome 10 (AD, missense) / 2 (AR, biallelic) / schwannomatosis (AD, LoF)",
    "CEP57": "MVA type 2 - the direct alternative to BUB1B",
    "TRIP13": "MVA type 3 - the direct alternative to BUB1B",
    "BUB1": "primary microcephaly 30; named in some MVA summaries",
    "BUB1B": "MVA type 1 - our hypothesis, printed for comparison not for credit",
}

#: Optional JSON object {gene: description} of further competing-diagnosis genes. It lives in
#: the quarantine and is never listed in this repository.
EXTRA_DIFFERENTIAL_FILE = REPO / "data/private/extra_differential.json"


def load_differential(extra: Path = EXTRA_DIFFERENTIAL_FILE) -> dict[str, str]:
    """Return the genes to check: the private extra genes first, then DIFFERENTIAL.

    The private genes come first so that the first unnamed gene in the output sits where it
    always has, and a rerun puts the summary of unnamed genes in the same place. A gene in
    both keeps its description from DIFFERENTIAL.

    Args:
        extra: JSON file holding an object of gene symbol -> description. An absent file
            adds nothing.

    Returns:
        The genes to check, in output order.
    """
    genes: dict[str, str] = {}
    if extra.is_file():
        import json

        loaded = json.loads(extra.read_text(encoding="utf-8"))
        if not isinstance(loaded, dict):
            raise ValueError(f"{extra.name} must hold a JSON object of gene symbol to description")
        genes.update({str(k): str(v) for k, v in loaded.items()})
    genes.update(DIFFERENTIAL)
    return genes


def no_extra_genes_note(genes: dict[str, str]) -> str | None:
    """Return the note to print when no gene beyond the shipped list was loaded, else None."""
    if len(genes) > len(DIFFERENTIAL):
        return None
    return (
        f"NOTE: no extra differential genes were loaded (data/private/extra_differential.json is absent "
        f"or adds none); only the {len(genes)} genes in the script were checked, so the count of "
        f"unnamed genes is absent."
    )


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

#: 🔴 Genes whose variants may appear at coordinate level in a PUBLIC artefact,
#: because they are the findings this submission actually reports and the
#: competition requires naming them. Everything else about this child stays in
#: the quarantine: a full run prints the proband's carrier status for unrelated
#: recessive conditions, which is individual-level health information that the
#: task does not need and a public repository must not carry. The privacy rules are in the README.
PUBLIC_SAFE_GENES = frozenset({"BUB1B", "LZTR1"})
REDACTED = "[redacted - not a reported finding; see the privacy rules in the README]"


def most_severe(bcsq: str) -> tuple[str | None, str | None]:
    best = gene = None
    for entry in bcsq.split(","):
        parts = entry.split("|")
        if len(parts) < 2:
            continue
        if "NMD_transcript" in parts[0] and best is not None:
            continue
        for candidate in SEVERITY:
            if candidate in parts[0]:
                if best is None or SEVERITY.index(candidate) < SEVERITY.index(best):
                    best, gene = candidate, parts[1]
                break
    return best, gene


def run(cmd: list[str]) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=False).stdout


def genome_wide_clinvar(
    vcf: Path, clinvar: Path, work: Path, min_dp: int, min_gq: int, min_mq: float, min_ab: float
):
    """Every high-quality ClinVar P/LP variant the proband actually carries."""
    plp: dict[tuple, tuple[str, str, str]] = {}
    regions = work / "plp_regions.txt"
    seen: set[tuple[str, str]] = set()
    with regions.open("w") as handle:
        rows = run(
            [
                "bcftools",
                "query",
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/GENEINFO\t%INFO/CLNDN\n",
                str(clinvar),
            ]
        )
        ordered = []
        for line in rows.splitlines():
            f = line.split("\t")
            if len(f) < 7 or "athogenic" not in f[4] or "Conflicting" in f[4]:
                continue
            plp[(f[0], f[1], f[2], f[3])] = (f[4], f[5], f[6])
            if (f[0], f[1]) not in seen:
                seen.add((f[0], f[1]))
                ordered.append((f[0], int(f[1])))
        for chrom, pos in sorted(ordered, key=lambda x: (x[0], x[1])):
            handle.write(f"{chrom}\t{pos}\n")
    print(f"ClinVar Pathogenic/Likely-pathogenic records: {len(plp):,}", file=sys.stderr)

    carried = []
    dropped_ab: list[tuple[str, str, str, str, float]] = []
    out = run(
        [
            "bcftools",
            "query",
            "-R",
            str(regions),
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\t[%DP]\t[%GQ]\t%INFO/MQ\t[%AD]\n",
            str(vcf),
        ]
    )
    for line in out.splitlines():
        f = line.split("\t")
        key = (f[0], f[1], f[2], f[3])
        if key not in plp or f[4] != "PASS":
            continue
        try:
            if int(f[6]) < min_dp or int(f[7]) < min_gq or float(f[8]) < min_mq:
                continue
            # Allele balance. DP/GQ/MQ alone do not catch a false heterozygote:
            # a call can sit at DP 34, GQ 99, MQ 58 and still be an artefact if
            # only 15% of its reads carry the alternate allele. Without this
            # filter an artefactual call in an unrelated gene entered the
            # "high-quality P/LP" list. Re-aligned from the raw reads it sat in a
            # region that recombines somatically in T cells, a known source of
            # chimeric reads in blood-derived DNA.
            ref_c, alt_c = (int(x) for x in f[9].split(",")[:2])
            total = ref_c + alt_c
            if total == 0 or alt_c / total < min_ab:
                dropped_ab.append((f[0], f[1], f[2], f[3], alt_c / total if total else 0.0))
                continue
        except (ValueError, IndexError):
            continue
        genotype = f[5].replace("|", "/")
        if genotype in ("0/0", "./."):
            continue
        carried.append((f[0], int(f[1]), f[2], f[3], genotype, f[6], plp[key]))
    regions.unlink(missing_ok=True)
    if dropped_ab:
        # Report, never silently discard: a rejected P/LP call is exactly the
        # kind of thing a reviewer must be able to see and disagree with.
        print(
            f"allele-balance filter rejected {len(dropped_ab)} ClinVar P/LP call(s) "
            f"below AB {min_ab}: "
            + ", ".join(f"{c}:{pos} {r}>{a} AB={ab:.3f}" for c, pos, r, a, ab in dropped_ab),
            file=sys.stderr,
        )
    return carried


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", type=Path, default=REPO / "data/private/raw/proband.vcf.gz")
    ap.add_argument("--clinvar", type=Path, default=REPO / "data/reference/clinvar/clinvar.vcf.gz")
    ap.add_argument("--tsv", type=Path, default=REPO / "data/private/work/damaging_genomewide.tsv")
    ap.add_argument("--min-dp", type=int, default=15)
    ap.add_argument("--min-gq", type=int, default=30)
    ap.add_argument("--min-mq", type=float, default=55.0)
    ap.add_argument(
        "--min-ab",
        type=float,
        default=0.25,
        help="minimum alternate-allele fraction for a heterozygous call. A true germline "
        "het sits near 0.5; below ~0.25 a call is far more likely an alignment artefact "
        "than a real allele. Default 0.25.",
    )
    ap.add_argument(
        "--public",
        action="store_true",
        help="🔴 Redact coordinates for every gene that is not a reported finding. "
        "Use this for any output committed to the repository or shown to anyone: a full "
        "run discloses the proband's carrier status for unrelated recessive conditions.",
    )
    args = ap.parse_args()

    if args.public:
        print(
            "# PUBLIC MODE - coordinates shown only for reported findings "
            f"({', '.join(sorted(PUBLIC_SAFE_GENES))}); all other loci redacted.\n"
        )

    work = REPO / "data" / "private" / "work"
    work.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("PART 1 - genome-wide ClinVar Pathogenic/Likely-pathogenic burden (no gene panel)")
    print("=" * 78)
    carried = genome_wide_clinvar(
        args.vcf, args.clinvar, work, args.min_dp, args.min_gq, args.min_mq, args.min_ab
    )
    by_gene: dict[str, list] = collections.defaultdict(list)
    for hit in carried:
        gene = hit[6][1].split(":")[0] if hit[6][1] and hit[6][1] != "." else "?"
        by_gene[gene].append(hit)
    print(f"\nThe proband carries {len(carried)} high-quality ClinVar P/LP variants genome-wide.\n")

    complete = {g: vs for g, vs in by_gene.items() if any(v[4] == "1/1" for v in vs) or len(vs) >= 2}
    print("Genes with a COMPLETE recessive P/LP genotype (homozygous, or two P/LP alleles):")
    if complete:
        for gene, vs in sorted(complete.items()):
            print(f"  {gene}: {len(vs)} allele(s)")
            for v in vs:
                if args.public and gene not in PUBLIC_SAFE_GENES:
                    print(f"     {REDACTED}")
                    continue
                print(f"     {v[0]}:{v[1]} {v[2]}>{v[3]} GT={v[4]} {v[6][0]} | {v[6][2][:60]}")
    else:
        print("  NONE. No second known-pathogenic allele anywhere in the genome.")
        print("  => no gene is completed by ClinVar-listed pathogenic variants alone.")

    print("\nHeterozygous single P/LP alleles (single alleles, or dominant findings to assess):")
    hidden = 0
    for gene, vs in sorted(by_gene.items()):
        if gene in complete:
            continue
        v = vs[0]
        if args.public and gene not in PUBLIC_SAFE_GENES:
            hidden += 1
            continue
        print(f"  {gene:<10} {v[0]}:{v[1]} {v[2]}>{v[3]} GT={v[4]} DP={v[5]}  {v[6][2][:58]}")
    if hidden:
        print(f"  + {hidden} further heterozygous P/LP allele(s) in genes unrelated to this")
        print("    phenotype. Deliberately NOT")
        print("    listed: they are individual-level health information about a child, the")
        print("    task does not require them, and this file is public. They are in the")
        print("    quarantined full run, and nowhere else.")

    print("\n" + "=" * 78)
    print("PART 2 - the phenotype-driven differentials, with evidence")
    print("=" * 78)

    differential = load_differential()
    note = no_extra_genes_note(differential)
    if note:
        # Visible in the output itself: run_all.sh discards standard error.
        print(note)
    per_gene: dict[str, list] = collections.defaultdict(list)
    for line in args.tsv.open():
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        consequence, gene = most_severe(f[9])
        if gene in differential:
            per_gene[gene].append((f[0], int(f[1]), f[2], f[3], consequence, f[5], f[6], f[8]))

    keys = [(c, p) for vs in per_gene.values() for c, p, *_ in vs]
    meta: dict[tuple, tuple[str, str]] = {}
    clnsig: dict[tuple, str] = {}
    if keys:
        regions = work / "diff_regions.txt"
        regions.write_text("\n".join(f"{c}\t{p}" for c, p in sorted(set(keys))) + "\n")
        for line in run(
            ["bcftools", "query", "-R", str(regions), "-f", "%CHROM\t%POS\t%REF\t%ALT\t%ID\n", str(args.vcf)]
        ).splitlines():
            f = line.split("\t")
            meta[(f[0], int(f[1]), f[2], f[3])] = (f[4], "")
        for line in run(
            [
                "bcftools",
                "query",
                "-R",
                str(regions),
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\n",
                str(args.clinvar),
            ]
        ).splitlines():
            f = line.split("\t")
            clnsig[(f[0], int(f[1]), f[2], f[3])] = f[4]
        regions.unlink(missing_ok=True)

    for gene, description in differential.items():
        variants = per_gene.get(gene, [])
        print(f"\n{gene} — {description}")
        if not variants:
            print("   no damaging variants at all  →  EXCLUDED at the genotype level")
            continue
        benign = 0
        for c, p, r, a, consequence, gt, dp, _gq in sorted(variants, key=lambda x: x[1]):
            key = (c, p, r, a)
            rsid = meta.get(key, (".", ""))[0]
            sig = clnsig.get(key, "-")
            if "enign" in sig:
                benign += 1
            novel = "novel" if rsid == "." else rsid
            if args.public and gene not in PUBLIC_SAFE_GENES:
                # Gene-level verdict is all a reader needs to check the exclusion.
                print(f"   {consequence:<14} GT={gt} {novel:<14} ClinVar={sig[:38]}   [coords redacted]")
                continue
            print(f"   {c}:{p} {r}>{a} {consequence:<14} GT={gt} DP={dp} {novel:<14} ClinVar={sig[:38]}")
        if benign == len(variants):
            print(f"   all {benign} are ClinVar-benign common polymorphisms  →  EXCLUDED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
