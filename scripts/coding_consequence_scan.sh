#!/usr/bin/env bash
#
# Genome-wide coding-consequence scan, one chromosome at a time.
#
# WHY THIS EXISTS
# ---------------
# Consequence annotation was the blocking gap for Track 1.
# The assumed route was a local Ensembl VEP install plus its ~25 GB GRCh38
# cache, which did not fit the machine we had. `bcftools csq` does the same job
# from an Ensembl GFF3 plus a reference FASTA, and it is haplotype-aware.
#
# The trick that makes it fit in a small box: never hold more than ONE
# chromosome's reference FASTA on disk at a time. Fetch, annotate, reduce to a
# compact TSV of damaging variants only, delete, move to the next chromosome.
#
#   peak extra disk  : one chromosome FASTA  (~250 MB, chr1)
#   peak RAM         : streaming only - bcftools pipes, nothing buffered
#   wall clock       : ~11 s of CPU per chromosome plus the download
#
# 🔴 Reads the quarantined VCF and writes only into data/private/. Nothing
# leaves the machine: the reference files are public downloads, and the proband
# VCF is never sent anywhere.
#
# Usage:  scripts/coding_consequence_scan.sh [output.tsv]
#
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="$REPO/data/reference/ensembl116"
GFF="$REF/Homo_sapiens.GRCh38.116.gff3.gz"
VCF="${MVA_VCF:-$REPO/data/private/raw/proband.vcf.gz}"
OUT="${1:-$REPO/data/private/work/damaging_genomewide.tsv}"
WORK="$(mktemp -d)"
BASE="https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna"

trap 'rm -rf "$WORK"' EXIT

for required in "$GFF" "$VCF"; do
  [ -s "$required" ] || { echo "missing: $required (run 'make refs' first)" >&2; exit 1; }
done
command -v bcftools >/dev/null || { echo "bcftools not on PATH" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
MANIFEST="${OUT%.tsv}.manifest.tsv"
: > "$OUT"
: > "$MANIFEST"

# Consequence classes we keep. INFO/ is required: bcftools csq defines BCSQ in
# both INFO and FORMAT, and an unqualified name is rejected as ambiguous.
FILT='INFO/BCSQ~"stop_gained" || INFO/BCSQ~"frameshift" || INFO/BCSQ~"splice_acceptor"
   || INFO/BCSQ~"splice_donor" || INFO/BCSQ~"start_lost" || INFO/BCSQ~"stop_lost"
   || INFO/BCSQ~"missense" || INFO/BCSQ~"inframe"'
FILT="${FILT//$'\n'/ }"

for C in {1..22} X Y; do
  FA="$WORK/c$C.fa"
  # 🔴 A skipped chromosome and a chromosome with no damaging variants both
  # produce no rows, and the output TSV cannot tell them apart. Downloads DO
  # fail transiently, so this used to `continue` with only a note on stderr --
  # which turns a transport error into a silent, permanent hole in a
  # "genome-wide" scan. Fail the run instead, and record what was processed in
  # a manifest next to the output so the claim is auditable after the fact.
  curl -sSf --max-time 900 -o "$FA.gz" "$BASE/Homo_sapiens.GRCh38.dna.chromosome.$C.fa.gz" \
    || { echo "chr$C: reference download FAILED - aborting rather than leaving a gap" >&2; exit 1; }
  gunzip -f "$FA.gz" || { echo "chr$C: gunzip FAILED - aborting rather than leaving a gap" >&2; exit 1; }
  bcftools faidx "$FA" 2>/dev/null || samtools faidx "$FA"

  # --ncsq 64: the default cap of 31 consequences per record silently truncates
  # genes with many annotated transcripts - BUB1B has 22.
  n=$(bcftools view -f PASS -r "$C" -Ou "$VCF" \
      | bcftools csq --ncsq 64 -l -f "$FA" -g "$GFF" -Ou - 2>"$WORK/err.$C" \
      | bcftools query -i "$FILT" \
          -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\t[%DP]\t[%AD]\t[%GQ]\t%INFO/BCSQ\n' - \
      | tee -a "$OUT" | wc -l)
  printf 'chr%-3s damaging=%-6s disk_avail=%s\n' "$C" "$n" "$(df -h / | awk 'NR==2{print $4}')"
  printf '%s\t%s\n' "$C" "$n" >> "$MANIFEST"
  rm -f "$FA" "$FA.fai" "$FA.gz"
done

printf 'contigs_processed\t%s\n' "$(wc -l < "$MANIFEST")" >> "$MANIFEST"
echo "SCAN_COMPLETE total=$(wc -l < "$OUT") contigs=$(grep -c . "$MANIFEST") -> $OUT"
echo "manifest -> $MANIFEST"

# 🔴 Note on FILTER: this scan reads `-f PASS` only. Non-PASS calls are
# deliberately excluded, so a damaging variant that fails a caller filter will
# not appear here. That is the right default (MQ40 and LowQual calls are mostly
# artefacts) but it is a CHOICE, not a completeness guarantee: chrY's only two
# damaging variants are `MQ40` and are invisible to this table by design.
# Interrogate a locus at every FILTER status separately before concluding a gene
# is variant-free.
