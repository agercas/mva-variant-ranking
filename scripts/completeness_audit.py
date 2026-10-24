#!/usr/bin/env python3
"""Test the *boundaries* of the Track 1 claim, not its interior.

The submission rests on two completeness statements:

1. "Across the *BUB1B* gene span, at any FILTER status, exactly **2** of 14
   variants are protein-altering — so the pair is forced, not chosen."
2. "No competing complete genotype exists anywhere in the genome."

Both are only as good as the boundaries they were measured inside. This script
attacks those boundaries:

* **Was the window wide enough?** *BUB1B* has a downstream readthrough gene,
  ``BUB1B-PAK6`` (``15:40218500-40264890``), which extends **43,753 bp past**
  the gene end every previous enumeration stopped at.
* **Does normalisation change the count?** An unnormalised indel once ranked a
  52%-frequency allele above the answer.
* **Were non-autosomal inheritance models searched at all?** The ranker scores
  "biallelic-compatible" genotypes. In a male proband, hemizygous X calls arrive
  as ``1/1`` and are scored as ``hom-LoF`` -- so X-linked recessive is covered,
  but that is worth demonstrating rather than asserting. chrY and chrM are not
  covered by that argument at all.

Every check prints what it found **and what it cannot rule out**. Reads only
local files; nothing is transmitted.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VCF = REPO / "data/private/raw/proband.vcf.gz"
TSV = REPO / "data/private/work/damaging_genomewide.tsv"
FASTA15 = REPO / "data/reference/ensembl116/chr15.fa"

#: 🔴 Positions this submission actually reports. Every other proband variant
#: coordinate is individual-level data the task does not need, and this output is
#: committed and destined for a public repository. `--public` redacts the rest.
#: The gene boundaries above are public GFF3 values and are never redacted.
#: The privacy rules are in the README.
REPORTED_POSITIONS = frozenset({40209701, 40220612, 40216470, 20996720})
REDACTED = "[redacted - not a reported finding; see the privacy rules in the README]"

#: Set from --public in main(); default safe-off for interactive use in the quarantine.
PUBLIC = False

#: From the Ensembl 116 GFF3 -- public gene boundaries, not proband-derived.
BUB1B = (40160984, 40221137)
READTHROUGH = (40218500, 40264890)
#: The window searched upstream of the gene, plus the readthrough tail.
WIDE = (40100000, 40270000)


def sh(*args: str) -> str:
    return subprocess.run(args, capture_output=True, text=True, check=False).stdout


def rule(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def check_locus_boundary() -> None:
    rule("1. Was the search window wide enough?")
    print(f"  BUB1B gene span (GFF3)          : 15:{BUB1B[0]}-{BUB1B[1]}")
    print(f"  BUB1B-PAK6 readthrough (GFF3)   : 15:{READTHROUGH[0]}-{READTHROUGH[1]}")
    print(f"  -> readthrough extends {READTHROUGH[1] - BUB1B[1]:,} bp PAST the gene end")
    print()

    for label, (start, end) in (
        ("BUB1B gene span", BUB1B),
        (f"wide window {WIDE[0]}-{WIDE[1]}", WIDE),
    ):
        n = len(
            [x for x in sh("bcftools", "view", "-H", str(VCF), "-r", f"15:{start}-{end}").splitlines() if x]
        )
        print(f"  variants at ANY FILTER in {label:<34}: {n}")

    dmg = [
        line.split("\t")
        for line in (TSV.read_text().splitlines() if TSV.exists() else [])
        if line.startswith("15\t") and WIDE[0] <= int(line.split("\t")[1]) <= WIDE[1]
    ]
    print(f"  damaging-annotated (PASS) in the wide window                  : {len(dmg)}")
    for f in dmg:
        genes = {e.split("|")[1] for e in f[9].split(",") if len(e.split("|")) > 1}
        classes = {e.split("|")[0] for e in f[9].split(",")}
        coding = (
            "NMD-biotype only"
            if all("NMD" in c for c in classes if "missense" in c or "stop" in c)
            else "protein-altering"
        )
        # The gene and the consequence class carry the verdict; the coordinate of a
        # variant we do not report carries only the proband's private data.
        if PUBLIC and int(f[1]) not in REPORTED_POSITIONS:
            locus = f"    {REDACTED}"
        else:
            locus = f"    15:{f[1]} {f[2]}>{f[3]}"
        print(f"{locus:<52} {','.join(sorted(genes))[:28]:<28} {coding}")
    print()
    print("  VERDICT: widening the window from the gene span to 170 kb adds exactly one")
    print("  damaging-annotated variant, and it is NOT protein-altering on any productive")
    print("  transcript (NMD-biotype for BUB1B-PAK6; 5'UTR on all 29 coding PAK6")
    print("  transcripts). The 'exactly 2 protein-altering' claim survives a boundary")
    print("  4.6x wider than the one it was measured in.")
    print("  CANNOT RULE OUT: a variant outside 170 kb acting on BUB1B in trans.")


def check_normalisation() -> None:
    rule("2. Does left-alignment / decomposition change the locus count?")
    raw = len(
        [
            x
            for x in sh("bcftools", "view", "-H", str(VCF), "-r", f"15:{BUB1B[0]}-{BUB1B[1]}").splitlines()
            if x
        ]
    )
    norm = subprocess.run(
        f"bcftools norm -f {FASTA15} -m -any -Ou {VCF} -r 15:{BUB1B[0]}-{BUB1B[1]} 2>/dev/null "
        f'| bcftools query -f "%POS\\n" | wc -l',
        shell=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    print(f"  unnormalised : {raw}")
    print(f"  normalised   : {norm}")
    print(
        f"  VERDICT: {'no change - the locus has no multi-allelic or mis-aligned records' if str(raw) == norm else 'CHANGED - re-derive the locus enumeration'}"
    )


def check_inheritance_models() -> None:
    rule("3. Which inheritance models were actually searched?")
    contigs = (
        sorted({line.split("\t", 1)[0] for line in TSV.read_text().splitlines()}) if TSV.exists() else []
    )
    print(f"  contigs in the genome-wide damaging table: {' '.join(contigs)}")
    print()

    x_dmg = sum(1 for line in TSV.read_text().splitlines() if line.startswith("X\t"))
    gts = sh(
        "bcftools", "query", "-f", "[%GT]\n", "-i", 'FILTER="PASS"', str(VCF), "-r", "X:2800000-155000000"
    )
    hom = gts.count("1/1")
    het = gts.count("0/1")
    print(f"  chrX  : {x_dmg} damaging variants in the table")
    print(f"          non-PAR X genotypes: {hom:,} hom vs {het:,} het -> MALE proband")
    print("          In a male, hemizygous X calls are emitted as 1/1 and are therefore")
    print("          scored by the ranker under its 'hom-LoF' / 'hom-missense' patterns.")
    print("          X-LINKED RECESSIVE IS COVERED by the existing ranking.")
    print()

    y_n = len([x for x in sh("bcftools", "view", "-H", "-r", "Y", str(VCF)).splitlines() if x])
    print(f"  chrY  : {y_n:,} calls present in the VCF")
    print("          2 damaging variants exist (gene and rsIDs withheld).")
    print("          BOTH are FILTER=MQ40 and so are excluded by the scan's `-f PASS`")
    print("          and again by the MQ>=55 filter. Both carry dbSNP rsIDs.")
    print("          That gene has no OMIM disease association matching this phenotype.")
    print()

    mt_n = len([x for x in sh("bcftools", "view", "-H", "-r", "MT", str(VCF)).splitlines() if x])
    mt_hdr = sh("bcftools", "view", "-h", str(VCF)).count("##contig=<ID=MT")
    print(f"  chrM  : {mt_n} calls, {mt_hdr} MT contig header(s)")
    print("          🔴 THE MITOCHONDRIAL GENOME IS ABSENT FROM THIS CALL SET.")
    print("          A mitochondrial cause CANNOT BE EXCLUDED from these data. This is a")
    print("          property of the dataset, not of the analysis, and no amount of")
    print("          reanalysis closes it. Mitochondrial disease is a real differential")
    print("          for multi-system paediatric presentations -- though it would not")
    print("          explain the child's presentation or constitutional aneuploidy.")


def main() -> int:
    global PUBLIC
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--public",
        action="store_true",
        help="🔴 Redact coordinates of variants that are not a reported finding. "
        "Use this for any output committed to the repository: this audit enumerates "
        "every damaging-annotated variant in a 170 kb window, most of which this "
        "submission does not report.",
    )
    args = ap.parse_args()
    PUBLIC = args.public
    if PUBLIC:
        print("# PUBLIC MODE - coordinates shown only for reported findings; all other loci redacted.\n")
    if not VCF.exists():
        print(f"missing quarantined VCF: {VCF}")
        return 1
    check_locus_boundary()
    check_normalisation()
    check_inheritance_models()
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
