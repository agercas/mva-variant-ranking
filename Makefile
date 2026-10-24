# Makefile for the public repository.
#
# Every target here is one step of the analysis; the README says which results file each
# one produces. Run `make help` for the list.
#
# Settings you may want to change:
#   PY=python3.12        which interpreter to use
#   VCF=<path>           the child's VCF (default: data/private/raw/*.vcf.gz)

VENV    := .venv
PY      ?= $(shell test -x $(VENV)/bin/python && echo $(VENV)/bin/python || echo python3)
RUFF    ?= $(shell test -x $(VENV)/bin/ruff && echo $(VENV)/bin/ruff || echo ruff)
SCRIPTS := scripts

.DEFAULT_GOAL := help

.PHONY: help
help: ## Show this help
	@grep -hE '^[a-zA-Z0-9_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: setup
setup: ## Install the package and the development dependencies into .venv
	python3 -m venv $(VENV)
	$(VENV)/bin/python -m pip install --upgrade pip
	$(VENV)/bin/python -m pip install -e ".[dev]"

# ------------------------------------------------------------ the analysis ---
REFDIR := data/reference
ENSEMBL_RELEASE := 116

.PHONY: refs
refs: ## Fetch public references (Ensembl GFF3 + ClinVar + gnomAD chr15 exomes, ~7.7 GB)
	@mkdir -p $(REFDIR)/ensembl$(ENSEMBL_RELEASE) $(REFDIR)/clinvar $(REFDIR)/gnomad
	@test -s $(REFDIR)/ensembl$(ENSEMBL_RELEASE)/Homo_sapiens.GRCh38.$(ENSEMBL_RELEASE).gff3.gz || \
		curl -sSf --max-time 900 -o $(REFDIR)/ensembl$(ENSEMBL_RELEASE)/Homo_sapiens.GRCh38.$(ENSEMBL_RELEASE).gff3.gz \
		https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.$(ENSEMBL_RELEASE).gff3.gz
	@test -s $(REFDIR)/clinvar/clinvar.vcf.gz || \
		curl -sSf --max-time 900 -o $(REFDIR)/clinvar/clinvar.vcf.gz \
		https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
	@test -s $(REFDIR)/clinvar/clinvar.vcf.gz.tbi || \
		curl -sSf --max-time 300 -o $(REFDIR)/clinvar/clinvar.vcf.gz.tbi \
		https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
	@test -s $(REFDIR)/gnomad/gnomad.exomes.v4.1.sites.chr15.vcf.bgz || \
		curl -sSf --max-time 3600 -o $(REFDIR)/gnomad/gnomad.exomes.v4.1.sites.chr15.vcf.bgz \
		https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr15.vcf.bgz
	@test -s $(REFDIR)/gnomad/gnomad.exomes.v4.1.sites.chr15.vcf.bgz.tbi || \
		curl -sSf --max-time 600 -o $(REFDIR)/gnomad/gnomad.exomes.v4.1.sites.chr15.vcf.bgz.tbi \
		https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr15.vcf.bgz.tbi
	@# chr15 and chr22 FASTA: needed by the from-first-principles consequence check
	@# (BUB1B and LZTR1 respectively). The genome-wide scan fetches the rest one
	@# chromosome at a time and deletes each after use.
	@for C in 15 22; do \
		test -s $(REFDIR)/ensembl$(ENSEMBL_RELEASE)/chr$$C.fa || { \
			curl -sSf --max-time 900 -o $(REFDIR)/ensembl$(ENSEMBL_RELEASE)/chr$$C.fa.gz \
			https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$$C.fa.gz && \
			gunzip -f $(REFDIR)/ensembl$(ENSEMBL_RELEASE)/chr$$C.fa.gz; }; \
	done
	@echo "references ready in $(REFDIR)/ (git-ignored)"

.PHONY: scan
scan: ## Genome-wide coding-consequence scan, one chromosome at a time (needs: refs, VCF)
	$(SCRIPTS)/coding_consequence_scan.sh

.PHONY: gnomad-af
gnomad-af: ## Real gnomAD AF for every damaging variant, genome-wide (replaces the rsID proxy)
	$(SCRIPTS)/gnomad_af_genomewide.sh

.PHONY: rank
rank: ## Rank genes by rare biallelic damaging genotype (phenotype-blind)
	$(PY) $(SCRIPTS)/rank_biallelic_genes.py $(ARGS)

.PHONY: sensitivity
sensitivity: ## Re-rank at every combination of the three filter thresholds (108 runs)
	$(SCRIPTS)/sensitivity.sh

.PHONY: differential
differential: ## Exclude competing diagnoses + list secondary findings (add ARGS="--public" to redact)
	$(PY) $(SCRIPTS)/differential_and_secondary_scan.py $(ARGS)

.PHONY: expanded-differential
expanded-differential: ## 43-gene MVA / chromosome-breakage / Fanconi / cancer-predisposition panel
	$(PY) $(SCRIPTS)/expanded_differential.py

.PHONY: verify-consequence
verify-consequence: ## Re-derive a consequence by hand from reference + GFF3 (see script --help)
	$(PY) $(SCRIPTS)/verify_consequence_from_first_principles.py $(ARGS)

.PHONY: aneuploidy
aneuploidy: ## Screen for mosaic whole-chromosome aneuploidy from BAF and depth
	$(PY) $(SCRIPTS)/aneuploidy_screen.py

.PHONY: overdispersion
overdispersion: ## BAF-overdispersion test for VARIEGATED aneuploidy (sensitive where per-chromosome means are not)
	$(PY) $(SCRIPTS)/baf_overdispersion.py

.PHONY: roh
roh: ## Test whether the locus sits in a run of homozygosity
	$(PY) $(SCRIPTS)/roh_check.py 2>/dev/null || $(PY) $(SCRIPTS)/roh_check.py

REF38DIR := data/reference/grch38
AMDIR    := data/reference/alphamissense
INSILICO := data/reference/insilico

.PHONY: refs38
refs38: ## Fetch the GRCh38 analysis set + AlphaMissense (~1.5 GB) for read-level work
	@mkdir -p $(REF38DIR) $(AMDIR)
	@test -s $(REF38DIR)/GRCh38.fa || { \
		test -s $(REF38DIR)/GRCh38_no_alt_plus_hs38d1.fna.gz || \
		curl -sSf --max-time 10800 -o $(REF38DIR)/GRCh38_no_alt_plus_hs38d1.fna.gz \
		"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"; \
		zcat $(REF38DIR)/GRCh38_no_alt_plus_hs38d1.fna.gz > $(REF38DIR)/GRCh38.fa; \
		samtools faidx $(REF38DIR)/GRCh38.fa; }
	@test -s $(AMDIR)/AlphaMissense_hg38.tsv.gz || \
		curl -sSfL --max-time 5400 -o $(AMDIR)/AlphaMissense_hg38.tsv.gz \
		https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
	@echo "read-level references ready (git-ignored)"

.PHONY: realign
realign: ## Targeted re-alignment of one locus from raw FASTQ: make realign ARGS="15 40110000 40270000 bub1b"
	@test -n "$(ARGS)" || (echo 'usage: make realign ARGS="<contig> <start> <end> <label>"'; exit 2)
	$(SCRIPTS)/targeted_realign.sh $(ARGS)

mito-refs: ## Fetch the rCRS / GRCh38 mitochondrial reference (~5 MB)
	@mkdir -p $(REFDIR)/mito
	@test -s $(REFDIR)/mito/MT.fa || { \
		curl -sSfL --max-time 300 -o $(REFDIR)/mito/MT.fa.gz \
		"https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"; \
		gunzip -f $(REFDIR)/mito/MT.fa.gz; }
	@echo "MT reference ready: $$(grep -vc '>' $(REFDIR)/mito/MT.fa) lines"

# Streams the child's FASTQ files and stores nothing. Give the two read files:
#   make mito-extract R1=<read 1>.fastq.gz R2=<read 2>.fastq.gz
mito-extract: mito-refs ## Recover chrM by STREAMING the raw FASTQs (~30 min, stores nothing)
	@test -n "$(R1)" -a -n "$(R2)" || (echo 'usage: make mito-extract R1=<read1.fastq.gz> R2=<read2.fastq.gz>'; exit 2)
	@$(PY) scripts/mito_from_fastq.py --mt-fasta $(REFDIR)/mito/MT.fa --files $(R1) $(R2) $(ARGS)

mito-call: ## Call and interpret mitochondrial variants from the streamed pileup
	@$(PY) scripts/mito_call_variants.py --mt-fasta $(REFDIR)/mito/MT.fa $(ARGS)

insilico-refs: ## Fetch the ClinGen-calibrated predictors for the BUB1B locus (~1.3 GB transient, 2.3 MB kept)
	@mkdir -p $(INSILICO) $(AMDIR)
	@test -s $(AMDIR)/AlphaMissense_hg38.tsv.gz || \
		curl -sSfL --max-time 5400 -o $(AMDIR)/AlphaMissense_hg38.tsv.gz \
		https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
	@test -s $(INSILICO)/bub1b_alphamissense.tsv || \
		zcat $(AMDIR)/AlphaMissense_hg38.tsv.gz | \
		awk -F'\t' '$$1=="chr15" && $$2>=40160984 && $$2<=40221137' > $(INSILICO)/bub1b_alphamissense.tsv
	@test -x $(INSILICO)/bigBedToBed || { \
		curl -sSfL --max-time 900 -o $(INSILICO)/bigBedToBed \
		https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed; \
		chmod +x $(INSILICO)/bigBedToBed; }
	@for tr in Vest PolyPhen2 Sift; do \
		test -s $(INSILICO)/bub1b_$$tr.bed || \
		$(INSILICO)/bigBedToBed -chrom=chr15 -start=40160984 -end=40221137 \
			"https://hgdownload.soe.ucsc.edu/gbdb/hg38/dbNsfp/dbNsfp$$tr.bb" \
			$(INSILICO)/bub1b_$$tr.bed; done
	@test -s $(INSILICO)/bub1b_revel.csv || { \
		curl -sSfL --max-time 3000 -o $(INSILICO)/revel.zip \
		"https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip"; \
		unzip -p $(INSILICO)/revel.zip revel_with_transcript_ids | \
		awk -F',' 'NR>1 && $$1=="15" && $$3!="." && $$3+0>=40160984 && $$3+0<=40221137' \
			> $(INSILICO)/bub1b_revel.csv; \
		rm -f $(INSILICO)/revel.zip; }
	@test -s $(INSILICO)/sac_family_genes.tsv || \
		zcat $(REFDIR)/ensembl$(ENSEMBL_RELEASE)/Homo_sapiens.GRCh38.$(ENSEMBL_RELEASE).gff3.gz | awk -F'\t' \
		'$$3=="gene" && $$9~/biotype=protein_coding/ {match($$9,/Name=[^;]+/); if(RSTART>0){n=substr($$9,RSTART+5,RLENGTH-5); if(n ~ /^(BUB1B|BUB1|BUB3|MAD1L1|MAD2L1|MAD2L1BP|TTK|CENPE|CEP57|TRIP13|SMC5|ESPL1|ZW10|ZWILCH|KNTC1|CDC20|AURKB|PLK1|NDC80|SPC24|SPC25|NUF2|CENPF|ERCC6L)$$/) print $$1"\t"$$4"\t"$$5"\t"n}}' \
		| sort -k4,4 > $(INSILICO)/sac_family_genes.tsv
	@test -s $(INSILICO)/sac_family_alphamissense.tsv || \
		$(PY) scripts/slice_alphamissense.py \
			--table $(AMDIR)/AlphaMissense_hg38.tsv.gz \
			--genes $(INSILICO)/sac_family_genes.tsv \
			--out $(INSILICO)/sac_family_alphamissense.tsv
	@echo "in-silico references ready (git-ignored). 🔴 Bulk downloads, queried locally: no coordinate is transmitted."

family-calibration: ## Calibrate AlphaMissense on the spindle-checkpoint gene family (needs: insilico-refs)
	@$(PY) scripts/family_calibration.py \
		--scores $(INSILICO)/sac_family_alphamissense.tsv \
		--genes $(INSILICO)/sac_family_genes.tsv $(ARGS)

insilico: ## Score BUB1B alleles against the ClinGen SVI calibrated predictors (needs: insilico-refs)
	@$(PY) scripts/insilico_panel.py $(ARGS)

.PHONY: geometry
geometry: ## N1002 in the public AlphaFold model (needs: pip install biopython)
	$(PY) $(SCRIPTS)/n1002_geometry.py --download data/reference/alphafold

.PHONY: alphamissense
alphamissense: ## AlphaMissense in context: the score, the same-site ladder, the calibrated band (public table, ~1 MB)
	$(PY) $(SCRIPTS)/alphamissense_bub1b.py --download $(AMDIR)

acmg: ## Formal ACMG/AMP classification of both alleles, on the Tavtigian point scale
	$(PY) $(SCRIPTS)/acmg_classify.py

# ---------------------------------------------------------- verification ---
.PHONY: hpo-refs crosscheck phenotype
hpo-refs: ## Fetch HPO disease-phenotype annotations and the ontology (~65 MB)
	$(MAKE) -C verification hpo-refs

crosscheck: ## Candidate-gene route: disease -> known genes -> compare (needs the VCF and network)
	$(MAKE) -C verification crosscheck

phenotype: ## Reverse route: does MVA1 match the child's features? (needs data/private/case_terms.json)
	$(MAKE) -C verification phenotype

# ------------------------------------------------------ the submission file ---
.PHONY: submission-csv
submission-csv: ## Rebuild and validate the predictions CSV into submission/
	$(PY) $(SCRIPTS)/build_track1_submission.py --out submission/bub1b-compound-het.csv

# --------------------------------------------------------------- quality ---
.PHONY: privacy-check
privacy-check: ## Scan the tree for data leakage and secrets
	$(PY) $(SCRIPTS)/privacy_guard.py --tree

.PHONY: lint
lint: ## Lint and format-check
	$(RUFF) check src tests scripts verification
	$(RUFF) format --check src tests scripts verification

.PHONY: test
test: ## Run the offline tests
	$(PY) -m pytest -q -m "not live"

.PHONY: check
check: privacy-check lint test ## Privacy scan, lint and tests
	@echo "all checks passed"
