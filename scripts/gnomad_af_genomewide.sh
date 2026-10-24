#!/usr/bin/env bash
# Fetch a real gnomAD v4.1 allele frequency for every damaging variant, genome-wide.
#
# The repo currently substitutes "has no dbSNP rsID in the proband VCF" for
# rarity everywhere except chr15, and measured that proxy to be right only 35%
# of the time. The stated reason was disk: a genome-wide gnomAD set is large.
# This box has 1.7 TB, so the proxy can simply be replaced with the real thing.
#
# Never holds more than one chromosome at a time: download, extract AF at our
# positions only, delete. No proband coordinate leaves the machine - the sites
# files are public bulk downloads read locally.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TSV="$REPO/data/private/work/damaging_genomewide.tsv"
OUT="$REPO/data/private/work/gnomad_af_genomewide.tsv"
WORK="$REPO/data/private/work/gnomad_tmp"
BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes"
mkdir -p "$WORK"; : > "$OUT"

for C in {1..22} X Y; do
  # positions we actually need on this chromosome
  awk -F'\t' -v c="$C" '$1==c {print "chr"c"\t"$2}' "$TSV" | sort -u -k2,2n > "$WORK/pos.txt"
  n_want=$(wc -l < "$WORK/pos.txt")
  if [ "$n_want" -eq 0 ]; then echo "chr$C: no damaging variants, skipped"; continue; fi

  F="$WORK/g.chr$C.vcf.bgz"
  if [ "$C" = "15" ] && [ -s "$REPO/data/reference/gnomad/gnomad.exomes.v4.1.sites.chr15.vcf.bgz" ]; then
    F="$REPO/data/reference/gnomad/gnomad.exomes.v4.1.sites.chr15.vcf.bgz"
  else
    curl -sSf --max-time 7200 -o "$F"     "$BASE/gnomad.exomes.v4.1.sites.chr$C.vcf.bgz"     || { echo "chr$C: download failed" >&2; rm -f "$F"; continue; }
    curl -sSf --max-time 900  -o "$F.tbi" "$BASE/gnomad.exomes.v4.1.sites.chr$C.vcf.bgz.tbi" || { echo "chr$C: tbi failed" >&2; rm -f "$F" "$F.tbi"; continue; }
  fi

  n=$(bcftools query -R "$WORK/pos.txt" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/nhomalt\t%INFO/AC\t%INFO/AN\n' "$F" 2>/dev/null | tee -a "$OUT" | wc -l)
  echo "chr$C: wanted $n_want positions, matched $n gnomAD records"
  [ "$F" = "$WORK/g.chr$C.vcf.bgz" ] && rm -f "$F" "$F.tbi"
done
echo "GNOMAD_AF_COMPLETE lines=$(wc -l < "$OUT") -> $OUT"
