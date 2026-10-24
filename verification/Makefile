# Self-contained: this directory does not touch the repo's root Makefile.
REPO := $(shell cd .. && pwd)
HPO  := $(REPO)/data/reference/hpo
PY   := python3

.PHONY: help hpo-refs crosscheck phenotype all
help:
	@grep -hE '^[a-z-]+:.*?##' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-14s\033[0m %s\n",$$1,$$2}'

hpo-refs: ## Fetch HPO disease-phenotype annotations and the ontology (~65 MB)
	@mkdir -p $(HPO)
	@test -s $(HPO)/phenotype.hpoa || curl -sSfL --max-time 600 -o $(HPO)/phenotype.hpoa \
		https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa
	@test -s $(HPO)/hp.obo || curl -sSfL --max-time 600 -o $(HPO)/hp.obo \
		https://purl.obolibrary.org/obo/hp.obo
	@echo "HPO references ready (git-ignored)."

crosscheck: ## Candidate-gene route: disease -> known genes -> compare (needs the VCF)
	@cd .. && $(PY) verification/crosscheck.py --public

phenotype: hpo-refs ## Reverse route: does MVA1 actually match this child's features?
	@cd .. && $(PY) verification/phenotype_reverse_check.py

all: crosscheck phenotype ## Both checks, into results/
