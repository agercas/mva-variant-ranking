# Walkthrough: reproduce the analysis from a clean machine

Run times are for 4 vCPUs, 16 GB of memory and 30 GB of disk. Each step names the
file in `results/` that holds its output. Steps 0 to 5 are enough to reach the answer.
They take about an hour, most of it spent downloading public references. The later steps
are the checks.

## The path at a glance

![From 4.7 million variants to two: each step with its tool, what it does, why, and what is left](pipeline.svg)

Section 3.1 of the report tells the same story in four stages and lists the changes.

## 0. Set up (3 minutes)

```bash
git clone https://github.com/agercas/mva-variant-ranking
cd mva-variant-ranking
make setup                     # Python environment with the test and lint tools
```

You also need `bcftools` (1.19 or later), and for the raw-read checks `minimap2`,
`samtools` and `bbmap`.

## 1. Get the data (2 minutes for the files the main path needs)

Request access to the gated dataset on Hugging Face, accept its rules, then:

```bash
hf download SageBio/mva-hackathon-2026-data --repo-type dataset \
    --include "*.vcf.gz*" "*.docx" --local-dir data/private/raw
```

This download is 317 MB and holds the VCF and the phenotype document. You need the eight FASTQ files
(about 85 GB) only for steps 7 and 8. Everything lands in `data/private/`,
which is git-ignored. For the phenotype check (step 6), copy the child's HPO IDs from the phenotype document into `data/private/case_terms.json` as a JSON
list.

## 2. Public references (25 minutes, about 7.7 GB)

```bash
make refs                      # Ensembl 116 gene models, ClinVar, gnomAD chr15, two chromosome FASTAs
```

## 3. Effect of every variant on the protein (20 minutes)

```bash
make scan                      # scripts/coding_consequence_scan.sh
```

The script handles one chromosome's reference at a time. It fetches the reference, annotates with `bcftools csq`, keeps the damaging
variants, deletes the reference and moves to the next chromosome. Peak extra disk is about 250 MB. It writes a manifest
so that a skipped chromosome cannot look like a chromosome with no hits.

## 4. Rank the genes (4 minutes)

```bash
make rank                                  # biallelic-ranking.txt
make rank ARGS="--no-clinvar --top 5"      # the same ranking with ClinVar withheld
make gnomad-af                             # real gnomAD frequency for every contig (slow, optional)
make sensitivity                           # 108 combinations of the filter thresholds (needs gnomad-af)
```

The script takes no gene, disease or phenotype as input. *BUB1B* ranks first of 175 genes in the
`results/biallelic-ranking.txt` and 41st of 175 without ClinVar.

## 5. Confirm the two variants (5 minutes)

```bash
make verify-consequence ARGS="--transcript ENST00000287598 --contig 15 --pos 40209701 \
    --ref T --alt G --fasta data/reference/ensembl116/chr15.fa --expect-protein-length 1050"
make differential expanded-differential    # competing diagnoses; the 43-gene panel
make aneuploidy overdispersion roh         # whole-chromosome gains and losses; homozygosity
python3 scripts/completeness_audit.py --public
```

`verify-consequence` rebuilds the coding sequence from the reference and translates it with
a hard-coded genetic code. Also run it on `--pos 40216470 --ref A --alt G`. It must say
that this position is not in the coding sequence.

## 6. The two independent routes (7 minutes)

```bash
make crosscheck                            # from the disease: known MVA genes, then what the child carries
make hpo-refs && make phenotype            # from the child's features: the rank of MVA1
```

## 7. Read-level checks (needs the FASTQ files)

```bash
make refs38
make realign ARGS="15 40110000 40270000 bub1b"            # exon depth, insert size, an independent re-call
make mito-extract R1=<read1.fastq.gz> R2=<read2.fastq.gz> # streams the reads, stores nothing, about 30 minutes
make mito-call
```

## 8. Interpret the variants (15 minutes, about 1.3 GB transient)

```bash
make acmg                                  # ACMG/AMP points for both alleles; no data needed
make alphamissense                         # public table; no data needed
make geometry                              # N1002 in the public AlphaFold model; needs biopython
make insilico-refs && make insilico        # REVEL, VEST4, PolyPhen-2, SIFT with controls
make family-calibration                    # generic threshold against this gene family
```

Read the full output. The calibrated predictors disagree, so the PP3 criterion (the ACMG criterion for
computer predictions) earns no points.

## 9. Build the predictions file (instant)

```bash
make submission-csv                        # writes and validates submission/*.csv
```

The builder checks the file against the organizers' format. It refuses a pair split over two
rows. It also scores the file against three possible answers with a copy of the organizers' rules.
Do not edit the CSV by hand.

## Everything at once

```bash
./run_all.sh
```

## The two independent checks (`verification/`)

Both scripts here start from a different place than the main ranking and share no code
with it. They sit in a separate folder, so the main analysis cannot import
them and they cannot import it.

| Check | Starts from | Question |
|---|---|---|
| `crosscheck.py` | the disease | Takes the MVA genes named in public sources, locates them with a different gene-coordinate source and transcript model, and finds the one that carries a known pathogenic allele in the child. It tests whether the main ranking's gene and variants are reached again. |
| `phenotype_reverse_check.py` | the child's features | Ranks every disease annotated in the Human Phenotype Ontology (HPO) by how well it matches the child's features, and reports the rank of MVA1. It tests whether the named disease is compatible with the child. |

### Check 1: starting from the disease

```bash
python3 verification/crosscheck.py --public
```

1. Take the gene list from Wikipedia, MedlinePlus and OMIM (nine genes).
2. Find each gene's position by symbol, through HGNC and the Ensembl REST service.
3. Read the child's variants in each gene straight from the VCF.
4. Ask which gene holds a variant that ClinVar classes as pathogenic or likely pathogenic.
5. In that gene, map every variant onto the Ensembl coding sequence (fetched by transcript ID) and
   translate it with a hard-coded codon table.
6. Place the two coding changes in the context of the protein, using UniProt.

It shares the VCF and `bcftools` with the main analysis, because there is one child and one
obvious tool. The gene list, the coordinate source, the transcript model, the protein annotation
and the logic are all different. **It starts from the disease**, so a wrong guess of MVA would pass
unnoticed. Check 2 covers that gap, and the genome-wide ranking is independent of the guess.

`--public` blanks the coordinates of variants that are not part of the reported finding.

### Check 2: starting from the child's features

```bash
make hpo-refs            # HPO files, public
# put the child's HPO IDs in data/private/case_terms.json, e.g. ["HP:0000001", "HP:0000002"]
python3 verification/phenotype_reverse_check.py
```

The child's HPO IDs come from the phenotype document in the dataset. They are read from a
git-ignored file and are not stored in this repository. The labels are read from `hp.obo`.

HPO is a graph, so every feature connects to every disease through some ancestor, up to the root
term `All`. A raw count of "features explained" would say every feature is explained, which means
nothing. So the script grades each link by how specific the connecting term is.
The measure is information content: the log of the number of diseases divided by the number that carry the term. The grades are:

* **direct**: the term is annotated on MVA1 itself
* **specific**: a shared ancestor held by fewer than about 2% of diseases (information content 4 or more)
* **generic**: information content between 2 and 4
* **vacuous**: below 2, which says almost nothing

Then the script ranks every disease by the summed information content of the terms it shares with the
child. It reports where MVA1 and a short list of rivals land.

This check shows whether MVA1 is compatible with the child. It starts from MVA, as the first check does,
so MVA1 as the answer needs the genotype evidence.
