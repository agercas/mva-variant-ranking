#!/usr/bin/env bash
#
# Targeted re-alignment of one locus straight from the raw FASTQ.
#
# WHY THIS EXISTS
# ---------------
# The dataset ships no BAM/CRAM, so anything that needs read-level evidence -
# copy number, structural variants, allele balance, insert size - was
# previously out of reach.
# Aligning 85 GB of FASTQ end to end takes ~10 h and is wasted effort when the
# question is about one gene.
#
# The trick: pull only the read pairs that could belong to the locus, then
# align that handful against the WHOLE genome.
#
#   1. Bait = the locus taken from the REPEAT-MASKED reference, so repeat
#      k-mers (50% of a typical locus) are N and match nothing. Without this
#      step the filter drags in a large fraction of the genome's reads.
#   2. bbduk k=31 exact match, keeping a pair if EITHER mate hits. Reads that
#      span a breakpoint still match on their anchored half, so split-read and
#      discordant-pair evidence survives the filter.
#   3. Align the survivors to the full decoy-aware genome - NOT to the locus.
#      Aligning to a slice would make every paralogous read a false variant;
#      against the whole genome they go where they belong.
#
# 🔴 Reads the quarantined FASTQ and writes only into data/private/. Nothing
# leaves the machine.
#
# What it was used for (report, sections 4 and 6):
#   BUB1B: flat depth across all 23 coding exons -> no exon-level deletion;
#   both calls reproduce at an alternate fraction of 0.543 / 0.517 through an
#   independent aligner; insert size median 437 bp, maximum 800 bp -> phase
#   cannot be read for two alleles 10,911 bp apart.
#
# Usage:
#   scripts/targeted_realign.sh <contig> <start> <end> <label>
# e.g.
#   scripts/targeted_realign.sh 15 40110000 40270000 bub1b
#
# Needs: bbduk.sh (bbmap) + java, minimap2, samtools, and the references from
# `make refs38`.

set -euo pipefail

CONTIG="${1:?usage: targeted_realign.sh <contig> <start> <end> <label>}"
START="${2:?}"
END="${3:?}"
LABEL="${4:?}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
W="$REPO/data/private/work"
REF38="$REPO/data/reference/grch38/GRCh38.fa"
FQ="${MVA_FASTQ_PREFIX:-$REPO/data/private/raw/proband_R}"
LANES="${MVA_LANES:-L001 L002 L003 L004}"
THREADS="${MVA_THREADS:-$(nproc)}"
MASKED="$W/chr${CONTIG}_rm.fa"
ENSEMBL="https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna"

mkdir -p "$W"
for required in "$REF38"; do
  [ -s "$required" ] || { echo "missing: $required (run 'make refs38')" >&2; exit 1; }
done
for tool in bbduk.sh minimap2 samtools; do
  command -v "$tool" >/dev/null || { echo "$tool not on PATH" >&2; exit 1; }
done

# ---- 1. repeat-masked bait -------------------------------------------------
if [ ! -s "$MASKED" ]; then
  echo "==> fetching repeat-masked chr$CONTIG"
  curl -sSf --max-time 900 -o "$MASKED.gz" \
    "$ENSEMBL/Homo_sapiens.GRCh38.dna_rm.chromosome.$CONTIG.fa.gz"
  gunzip -f "$MASKED.gz"
  samtools faidx "$MASKED"
fi

BAIT="$W/bait_$LABEL.fa"
samtools faidx "$MASKED" "$CONTIG:$START-$END" > "$BAIT"
python3 - "$BAIT" <<'PY'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
n = seq.upper().count("N")
print(f"    bait {len(seq)} bp, {100*n/len(seq):.0f}% repeat-masked, {len(seq)-n} bp usable")
PY

# ---- 2. k-mer prefilter over the raw FASTQ ---------------------------------
: > "$W/${LABEL}_R1.fq"
: > "$W/${LABEL}_R2.fq"
for L in $LANES; do
  echo "==> baiting lane $L"
  bbduk.sh -Xmx24g "threads=$THREADS" \
    "in1=${FQ}_${L}_R1_001.fastq.gz" "in2=${FQ}_${L}_R2_001.fastq.gz" \
    "ref=$BAIT" k=31 mm=f rcomp=t \
    "outm=$W/_m_R1.fq" "outm2=$W/_m_R2.fq" ow=t 2>&1 | grep -E 'Contaminants:' || true
  cat "$W/_m_R1.fq" >> "$W/${LABEL}_R1.fq"
  cat "$W/_m_R2.fq" >> "$W/${LABEL}_R2.fq"
  rm -f "$W/_m_R1.fq" "$W/_m_R2.fq"
done
PAIRS=$(( $(wc -l < "$W/${LABEL}_R1.fq") / 4 ))
echo "    recovered $PAIRS read pairs"

# ---- 3. align the survivors to the WHOLE genome ----------------------------
echo "==> aligning to the full decoy-aware genome"
minimap2 -t "$THREADS" -ax sr -R "@RG\tID:prb\tSM:PROBAND\tPL:ILLUMINA" \
  "$REF38" "$W/${LABEL}_R1.fq" "$W/${LABEL}_R2.fq" 2>/dev/null \
  | samtools sort -@ 4 -o "$W/$LABEL.bam" -
samtools index "$W/$LABEL.bam"

echo
samtools flagstat "$W/$LABEL.bam" | head -5
echo
echo "TARGETED_REALIGN_COMPLETE pairs=$PAIRS -> $W/$LABEL.bam"
echo
echo "Interpreting the output:"
echo "  * Depth is ONLY meaningful over unique sequence. The bait was"
echo "    repeat-masked, so repeat windows are empty by construction, not by"
echo "    biology - check the repeat fraction of any interval before reading"
echo "    its coverage (all 23 BUB1B CDS exons are >80% unique; 2 kb windows"
echo "    across the locus mostly are not)."
echo "  * Soft-clip clusters at repeat boundaries are an artefact of the same"
echo "    masking. A real breakpoint sits in unique sequence."
