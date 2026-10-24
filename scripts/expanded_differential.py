"""Exclude competing diagnoses across a 43-gene panel.

The panel has four groups: genes reported to cause mosaic variegated aneuploidy
(MVA) or related chromosome-instability disorders, chromosome-breakage syndromes
(Bloom, Nijmegen and similar), the Fanconi anaemia pathway, and genes that
predispose to childhood cancers. The last three groups matter because severe
growth restriction, short stature and a childhood solid tumour are also the
classic presentation of the chromosome-breakage syndromes.

For each gene the script prints how many variants the child carries in the gene
span, how many have no dbSNP ID, and which ClinVar Pathogenic/Likely-pathogenic
alleles the child carries. Local files only; nothing is transmitted.
"""

from __future__ import annotations

import gzip
import re
import subprocess
from collections import defaultdict

VCF = "data/private/raw/proband.vcf.gz"
GFF = "data/reference/ensembl116/Homo_sapiens.GRCh38.116.gff3.gz"
CLINVAR = "data/reference/clinvar/clinvar.vcf.gz"

PANEL = {
    "MVA / aneuploidy": ["BUB1B", "CEP57", "TRIP13", "BUB1", "MAD2L1BP", "SMC5", "CENPE", "PLK4"],
    "chromosome breakage": ["BLM", "NBN", "ATM", "LIG4", "RECQL4", "ESCO2", "MCM4", "ATR", "RNU4ATAC"],
    "Fanconi": [
        "FANCA",
        "FANCB",
        "FANCC",
        "FANCD2",
        "FANCE",
        "FANCF",
        "FANCG",
        "FANCI",
        "FANCL",
        "FANCM",
        "BRCA1",
        "BRCA2",
        "BRIP1",
        "PALB2",
        "RAD51C",
        "SLX4",
        "ERCC4",
        "UBE2T",
    ],
    "RMS / cancer predisposition": ["TP53", "DICER1", "NF1", "RB1", "HRAS", "PTPN11", "SMARCB1", "SMARCA4"],
}
WANTED = {g for v in PANEL.values() for g in v}

# 1. gene spans from the GFF3
spans: dict[str, tuple[str, int, int]] = {}
with gzip.open(GFF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"Name=([^;]+)", f[8])
        if m and m.group(1) in WANTED and re.fullmatch(r"[0-9XYM]+", f[0]):
            spans[m.group(1)] = (f[0], int(f[3]), int(f[4]))

missing = WANTED - spans.keys()
if missing:
    print(f"[gene span not found in GFF3: {sorted(missing)}]\n")

# 2. ClinVar P/LP index for those spans only
clinvar: dict[tuple[str, int, str, str], str] = {}
for _gene, (c, s, e) in spans.items():
    out = subprocess.run(
        [
            "bcftools",
            "query",
            "-r",
            f"{c}:{s}-{e}",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\n",
            CLINVAR,
        ],
        capture_output=True,
        text=True,
    ).stdout
    for line in out.splitlines():
        ch, pos, ref, alt, sig = line.split("\t")
        if "athogenic" in sig and "Conflicting" not in sig:
            clinvar[(ch, int(pos), ref, alt)] = sig

print(f"ClinVar P/LP records inside the panel spans: {len(clinvar)}\n")

# 3. the proband's variants in each span
findings = defaultdict(list)
for gene, (c, s, e) in sorted(spans.items()):
    out = subprocess.run(
        [
            "bcftools",
            "query",
            "-r",
            f"{c}:{s}-{e}",
            "-f",
            "%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%DP]\t%ID\n",
            VCF,
        ],
        capture_output=True,
        text=True,
    ).stdout
    n_tot = n_novel = 0
    for line in out.splitlines():
        pos, ref, alt, filt, gt, dp, rsid = line.split("\t")
        n_tot += 1
        if rsid == ".":
            n_novel += 1
        hit = clinvar.get((c, int(pos), ref, alt))
        if hit:
            findings[gene].append((pos, ref, alt, gt, dp, filt, rsid, hit))
    flag = ""
    if findings[gene]:
        flag = f"  <-- {len(findings[gene])} ClinVar P/LP"
    print(f"{gene:<10} chr{c:<3} {n_tot:>5} variants  {n_novel:>3} novel{flag}")

print("\n=== ClinVar Pathogenic/Likely-pathogenic alleles carried, by gene ===")
if not findings:
    print("  none")
for gene, rows in sorted(findings.items()):
    print(f"\n{gene}:")
    for pos, ref, alt, gt, dp, filt, rsid, sig in rows:
        zyg = "HOM" if gt in ("1/1", "1|1") else "het"
        print(f"   {pos} {ref}>{alt} {zyg} DP={dp} {filt} {rsid} {sig}")
