#!/usr/bin/env bash
#
# Rerun every check and write the public view of each result into results/.
#
# Before you run it
#   make setup                  # Python environment
#   make refs                   # Ensembl GFF3, ClinVar, gnomAD chr15 and two chromosome FASTAs (~7.7 GB)
#   the gated VCF in data/private/raw/   (see docs/WALKTHROUGH.md, step 1)
#
# Optional, for the checks that need more disk or the raw reads. Each has its own
# make target and each writes its own file in results/:
#   make gnomad-af              # real gnomAD frequencies for every contig: adds the tighter ranking
#   make hpo-refs               # HPO files for the phenotype check
#   data/private/case_terms.json  # the child's HPO IDs, one JSON list (never published)
#   make insilico-refs && make insilico        -> insilico-panel.txt
#   make family-calibration                    -> family-calibration.txt
#   make mito-extract R1=... R2=... && make mito-call   -> mitochondrial-screen.txt (view: mito)
#
# Files that describe the child go through scripts/public_view.py, so other genes'
# genotypes, rsIDs, term names and mitochondrial positions never reach results/.
#
# This script does not upload anything.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$HERE"
while [ ! -f "$REPO/pyproject.toml" ] && [ "$REPO" != "/" ]; do REPO="$(dirname "$REPO")"; done
[ -f "$REPO/pyproject.toml" ] || { echo "cannot find the repository root above $HERE" >&2; exit 1; }

SUBMIT="$HERE/submission"

OUT="$HERE/results"
STAMP="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
WORK="data/private/work"
cd "$REPO"
mkdir -p "$OUT" "$SUBMIT"
echo "repository : $REPO"
echo "results    : $OUT"

step() { printf '\n==> %s\n' "$*"; }

# emit <file> <header> <command...>        output goes to results/ as it is
emit() {
  local out="$1" header="$2"; shift 2
  { echo "# $header"; echo "# Regenerated $STAMP by run_all.sh"; echo; "$@" 2>&1; } > "$OUT/$out"
}
# emit_public <file> <kind> <header> <command...>    output goes through public_view.py first
emit_public() {
  local out="$1" kind="$2" header="$3"; shift 3
  local raw; raw="$(mktemp)"
  { echo "# $header"; echo "# Regenerated $STAMP by run_all.sh"; echo; "$@" 2>/dev/null; } > "$raw"
  python3 scripts/public_view.py --kind "$kind" < "$raw" \
    | awk -v note="# Public view: this file went through scripts/public_view.py --kind $kind; other genes' genotypes, rsIDs, term names and mitochondrial positions are withheld." \
        'BEGIN{done=0} { if (!done && $0 !~ /^#/) { print note; done=1 } print }' > "$OUT/$out"
  rm -f "$raw"
}

step "1/10  coding consequence of every PASS variant, one chromosome at a time (~20 min)"
if [ -s "$WORK/damaging_genomewide.tsv" ]; then
  echo "    reusing $WORK/damaging_genomewide.tsv (delete it to rerun)"
else
  scripts/coding_consequence_scan.sh
fi

step "2/10  ranking genes by recessive-compatible genotype, and the ablation without ClinVar"
emit_public biallelic-ranking.txt ranking "No gene, no panel and no phenotype is an input to this ranking." \
  bash -c '
    echo "### GENOME-WIDE (real gnomAD AF on chr15; rsID rarity proxy elsewhere) ###"
    python3 scripts/rank_biallelic_genes.py --top 15
    echo; echo "### chr15 ONLY, under a REAL gnomAD allele-frequency filter ###"
    python3 scripts/rank_biallelic_genes.py --contig 15 --top 10
    echo; echo "### ABLATION: genotype pattern alone, ALL ClinVar evidence withheld ###"
    python3 scripts/rank_biallelic_genes.py --no-clinvar --top 5'
if [ -s "$WORK/gnomad_af_genomewide.tsv" ]; then
  emit_public biallelic-ranking-real-af.txt ranking "The same ranking with a real gnomAD v4.1 frequency for every contig." \
    bash -c '
      echo "### REAL gnomAD AF ON EVERY CONTIG ###"
      python3 scripts/rank_biallelic_genes.py --gnomad-tsv data/private/work/gnomad_af_genomewide.tsv --top 10
      echo; echo "### ABLATION: ClinVar withheld ###"
      python3 scripts/rank_biallelic_genes.py --gnomad-tsv data/private/work/gnomad_af_genomewide.tsv --no-clinvar --top 5'
  emit sensitivity-sweep.txt "BUB1B rank at every combination of the filter thresholds (108 runs)." \
    scripts/sensitivity.sh
else
  echo "    skipping the real-frequency ranking and the 108-run sweep: run 'make gnomad-af' first"
fi

step "3/10  consequences re-derived by hand from the reference and the GFF3"
emit consequence-verification.txt "Hard-coded genetic code. Does not trust bcftools csq." \
  bash -c '
    for spec in "ENST00000287598 15 40209701 T G 1050" "ENST00000287598 15 40220612 T G 1050" "ENST00000646124 22 20996720 C G 0"; do
      set -- $spec
      echo "-----------------------------------------------------------------"
      if [ "$6" = "0" ]; then
        python3 scripts/verify_consequence_from_first_principles.py --transcript "$1" --contig "$2" --pos "$3" --ref "$4" --alt "$5" --fasta data/reference/ensembl116/chr"$2".fa
      else
        python3 scripts/verify_consequence_from_first_principles.py --transcript "$1" --contig "$2" --pos "$3" --ref "$4" --alt "$5" --fasta data/reference/ensembl116/chr"$2".fa --expect-protein-length "$6"
      fi
      echo
    done
    echo "-----------------------------------------------------------------"
    echo "NEGATIVE CONTROL - the intronic variant must NOT be in the CDS:"
    python3 scripts/verify_consequence_from_first_principles.py --transcript ENST00000287598 --contig 15 --pos 40216470 --ref A --alt G --fasta data/reference/ensembl116/chr15.fa || true'

step "4/10  competing diagnoses: all of ClinVar, then the named differentials"
emit_public differential-exclusion.txt differential "Genome-wide ClinVar P/LP burden, no gene panel, then the named differentials." \
  python3 scripts/differential_and_secondary_scan.py --public

step "5/10  the 43-gene panel"
emit_public expanded-differential.txt expanded "MVA, chromosome-breakage, Fanconi and cancer-predisposition genes." \
  python3 scripts/expanded_differential.py

step "6/10  mosaic aneuploidy, its variance version, and runs of homozygosity"
emit aneuploidy-screen.txt "Per-chromosome allele balance and depth, with a detection limit." python3 scripts/aneuploidy_screen.py
emit overdispersion.txt "Variance test for mosaic aneuploidy." python3 scripts/baf_overdispersion.py
emit roh-check.txt "Is the gene in a run of homozygosity?" python3 scripts/roh_check.py

step "7/10  completeness audit: window, normalisation, inheritance models"
emit_public completeness-audit.txt completeness "Attacks the edges of the claim rather than its interior." \
  python3 scripts/completeness_audit.py --public

step "8/10  independent cross-checks"
emit crosscheck.txt "Disease first: known MVA genes, then what the child carries. Independent of the ranking." \
  python3 verification/crosscheck.py --public
if [ -s data/reference/hpo/hp.obo ] && [ -s data/private/case_terms.json ]; then
  emit_public phenotype-reverse-check.txt phenotype "Does MVA1 match the child's features?" \
    python3 verification/phenotype_reverse_check.py
else
  echo "    skipping the phenotype check: needs 'make hpo-refs' and data/private/case_terms.json"
fi

step "9/10  ACMG classification and AlphaMissense in context (no data needed)"
emit acmg-classification.txt "Both alleles on the ACMG/AMP point scale, and what would move allele 2." \
  python3 scripts/acmg_classify.py
emit alphamissense-bub1b.txt "AlphaMissense on BUB1B: the child's missense allele in context." \
  python3 scripts/alphamissense_bub1b.py --download data/reference/alphamissense
if python3 -c "import Bio" 2>/dev/null; then
  emit n1002-geometry.txt "Residue N1002 in the public AlphaFold model of BubR1." \
    python3 scripts/n1002_geometry.py --download data/reference/alphafold
else
  echo "    skipping the structure measurement: pip install biopython"
fi

step "10/10  rebuild and validate the predictions CSV, then scan the package"
python3 scripts/build_track1_submission.py --out "$SUBMIT/bub1b-compound-het.csv"
python3 scripts/privacy_guard.py --paths "$HERE"

printf '\ndone. Results are in %s\n' "$OUT"
