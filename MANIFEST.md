# Contents of this repository

This repository holds the code and results behind the report, and only what the answer needs. This file
lists every file in this repository.

## Structure

| Path | What it is |
|---|---|
| `README.md` | Start here: the answer, how to read the results and how to rerun |
| `MANIFEST.md` | This file |
| `LICENSE` | CC BY 4.0, the licence the competition rules require |
| `docs/pipeline.svg` | The diagram of the analysis: at each step, the tool, what it does and why |
| `docs/references.md` | Every tool, dataset and paper used, with links |
| `docs/WALKTHROUGH.md` | The steps, in order, with run times and the file each one writes |
| `results/` | The output of every check (see "Results" below) |
| `submission/` | The four files uploaded to the competition: predictions CSV, report, methods description (`.xlsx` and `.md`) |
| `run_all.sh` | Reruns every check on a machine that has the data |

## Code

### The answer (steps 1 to 5 of the report)

| File | What it does | Writes |
|---|---|---|
| `scripts/coding_consequence_scan.sh` | Works out the effect on the protein of every PASS variant with `bcftools csq`, one chromosome at a time (peak extra disk about 250 MB) | `damaging_genomewide.tsv` (kept private) |
| `scripts/gnomad_af_genomewide.sh` | Looks up the gnomAD v4.1 frequency of each damaging variant | (private table) |
| `scripts/rank_biallelic_genes.py` | Ranks genes by whether they carry a rare genotype that fits recessive inheritance. It takes no gene, disease or phenotype as input | `biallelic-ranking.txt` |
| `scripts/sensitivity.sh` | Repeats the ranking for 108 combinations of the filter thresholds | `sensitivity-sweep.txt` |
| `scripts/verify_consequence_from_first_principles.py` | Rebuilds the coding sequence from the reference and translates it with a hard-coded genetic code. It does not rely on the annotation tool | `consequence-verification.txt` |
| `scripts/completeness_audit.py` | Checks the edges of the search: window size, variant normalisation, and coverage of X, Y and mitochondria | `completeness-audit.txt` |
| `scripts/targeted_realign.sh` | Re-aligns the reads around one gene from the raw FASTQ files. It gives exon depth, insert size and an independent check of the calls | (private) |

### Excluding other explanations (report section 4)

| File | What it does | Writes |
|---|---|---|
| `scripts/differential_and_secondary_scan.py` | Finds which of the child's variants match any ClinVar pathogenic or likely pathogenic record, then checks named competing diagnoses | `differential-exclusion.txt` |
| `scripts/expanded_differential.py` | The same for a 43-gene panel (MVA, chromosome breakage, Fanconi pathway, cancer predisposition) | `expanded-differential.txt` |
| `scripts/aneuploidy_screen.py` | Looks for whole-chromosome gains or losses in allele balance and depth | `aneuploidy-screen.txt` |
| `scripts/baf_overdispersion.py` | A variance test for mosaic aneuploidy | `overdispersion.txt` |
| `scripts/roh_check.py` | Checks whether the gene sits in a run of homozygosity | `roh-check.txt` |
| `scripts/mito_from_fastq.py` | Rebuilds the mitochondrial genome by streaming the raw reads and storing nothing | (private pileup) |
| `scripts/mito_call_variants.py` | Calls mitochondrial variants that rise above the background noise level measured from the data | `mitochondrial-screen.txt` |
| `verification/crosscheck.py` | Starts from the disease, takes the genes named in public sources and checks what the child carries | `crosscheck.txt` |
| `verification/phenotype_reverse_check.py` | Ranks all diseases in HPO (the Human Phenotype Ontology) by how well they match the child's features, and reports where MVA1 falls | `phenotype-reverse-check.txt` |
| `verification/Makefile` | Targets for the two checks above | |

### Interpreting the two variants (report section 5)

| File | What it does | Writes |
|---|---|---|
| `scripts/acmg_classify.py` | Applies the ACMG/AMP criteria to both alleles and shows what would change the class of each | `acmg-classification.txt` |
| `scripts/alphamissense_bub1b.py` | Reads the public AlphaMissense table for BUB1B: the score, the scores of the other changes at the same site, and the ClinGen band | `alphamissense-bub1b.txt` |
| `scripts/n1002_geometry.py` | Measures residue N1002 in the public AlphaFold model: how buried it is, its hydrogen bonds and its distance from the lysine at the ATP site | `n1002-geometry.txt` |
| `scripts/insilico_panel.py` | Scores the missense allele with REVEL, VEST4, PolyPhen-2 and SIFT, with positive and negative controls | `insilico-panel.txt` |
| `scripts/family_calibration.py` | Tests whether the general predictor thresholds fit the spindle-checkpoint gene family | `family-calibration.txt` |
| `scripts/slice_alphamissense.py` | Cuts the AlphaMissense genome table down to the gene family (the input to the script above) | (public slice) |

### The predictions file and the scorer

| File | What it does |
|---|---|
| `scripts/build_track1_submission.py` | Writes and validates the predictions CSV, then scores it against three possible answers using a copy of the organizers' rules |
| `src/mva/scoring.py` | A copy of the organizers' scoring rules (rank points and F-max), checked against their `evaluation.py` |
| `src/mva/genomics/coordinates.py` | Variant coordinates; converts between the no-`chr` VCF names and the `chr` names the submission needs |
| `src/mva/genomics/submission.py` | The submission row format and its validator |
| `src/mva/genomics/__init__.py`, `src/mva/__init__.py` | Package files |
| `scripts/methods_xlsx_to_md.py` | Makes the Markdown methods description from the filled-in Excel template, so the two cannot differ |

### Checks and tests

| File | What it does |
|---|---|
| `scripts/privacy_guard.py` | Blocks genotype files, reads and credentials from entering the repository |
| `scripts/public_view.py` | Turns a raw result file into a version that is safe to publish. It withholds other genes' genotypes, rsIDs, term names and mitochondrial positions |
| `tests/conftest.py` | Shared test setup (synthetic variants only) |
| `tests/test_coordinates.py`, `tests/test_scoring.py`, `tests/test_submission.py` | Coordinate handling, the scorer, and the predictions file |
| `tests/test_privacy_guard.py` | The privacy guard |
| `tests/test_public_view.py` | The public-view tool, on synthetic text |
| `tests/test_verification_phenotype.py` | The reverse phenotype check, on a synthetic ontology |

## Results

Each file in `results/` is the output of one script above. Files that describe the
child go through `scripts/public_view.py` first. [`results/RECORDED-CHECKS.md`](results/RECORDED-CHECKS.md)
lists the numbers in the report that come from runs whose output files are not yet in
`results/`. For each number it gives the command that regenerates it.

## Also in the tree

| File | What it is |
|---|---|
| `.gitignore` | Keeps genotype files, read files and `data/` out of the repository |
| `Makefile` | One target per step of the analysis, plus `setup`, `refs` and `check`; run `make help` for the list |
| `pyproject.toml` | Package metadata, dependencies, and the pytest and ruff settings |

† Some result files named in `run_all.sh` (`sensitivity-sweep.txt`, `expanded-differential.txt`, `overdispersion.txt`, `roh-check.txt`, `biallelic-ranking-real-af.txt`) are written only when the pipeline is rerun on the data machine. Until then, their numbers are listed in `results/RECORDED-CHECKS.md`.

The scripts read the child's genotype file from `data/private/raw/proband.vcf.gz` (or the `MVA_VCF` variable). Place or link the file there.
