# Numbers in the report whose output files are not in this folder yet

Every number in the report comes from one of three places:

1. **A file in this folder.** Open it and read the number.
2. **A public table that anyone can fetch.** The script that reads it is in `scripts/`.
3. **A run listed below.** These runs were made on the machine that holds the child's data
   and written down in the authors' working notes. Their output files are not yet in
   `results/`. The command in the last column regenerates each one;
   a reader can check the numbers once that output is committed.

`run_all.sh` writes most of these files when it is run on a machine with the data. The
rest need the raw reads or extra downloads and have their own `make` targets.

| # | What the report says | Value | Regenerate with |
|---|---|---|---|
| R1 | Size of the call set | 4,740,790 PASS records | `bcftools view -H -f PASS <vcf> \| wc -l` |
| R2 | The gene ranking with a real gnomAD frequency on every contig | BUB1B first of **27** genes (the pool shrank 209 → 48 → 40 → 27 as the real frequencies, an annotation fix and an allele-balance filter were added; 209 is the first run with the dbSNP proxy, and the committed run with the fix and the filter but the proxy gives 175). The same run found that "no dbSNP ID" agrees with the real frequency in 26.8% of cases genome-wide (227 of 847), and that an annotation fix removed 1,593 of 17,558 damaging records (9.07%) that were never protein-changing | `make gnomad-af && python3 scripts/rank_biallelic_genes.py --gnomad-tsv data/private/work/gnomad_af_genomewide.tsv` |
| R3 | The ranking at 108 combinations of the filter thresholds | BUB1B first in all 108; candidate pool 14 to 64 genes | `make sensitivity` |
| R4 | The ranking with no filters | BUB1B at position 584; the top of the list is made of mapping artefacts (mucin, HLA and keratin-associated genes) | `python3 scripts/rank_biallelic_genes.py` with the three filters switched off |
| R5 | Splice effect of every variant in the locus | SpliceAI 1.3.1 (run on the local machine): all 15 variants at or below 0.03; the deep intronic variant in row 3 of the predictions scores 0.00; a synthetic splice-donor control scores 0.98 | `spliceai -I <locus vcf> -R <chr15.fa> -A grch38 -D 500` |
| R6 | Exon-level deletion | 22,234 read pairs recovered from the raw reads (99.0% mapped); all 23 coding exons at 0.77 to 1.44 times the median depth of 53.5×; no exon at half depth | `make refs38 && make realign ARGS="15 40110000 40270000 bub1b"` |
| R7 | Independent check of the two calls | Re-aligned with minimap2 and a fresh pileup: alternate-allele fractions 0.543 (stop-gain) and 0.517 (missense) | same as R6 |
| R8 | Read-pair distances at the locus | Insert size median 437 bp, 99th percentile 779 bp, maximum 800 bp. Every gap between neighbouring heterozygous sites at the locus is longer than 800 bp (smallest gap involving either variant 4,142 bp, largest 16,809 bp) | same as R6 |
| R9 | Upstream regulatory region | 15:40100000-40161000 holds 88 variants, 7 with no dbSNP ID, all 7 indels in repeats with the usual artefact signature; the site 44 kb upstream of the transcription start that was reported in Japanese families carries common SNPs only | `bcftools view -r 15:40100000-40161000 <vcf>` |
| R10 | Number of ClinVar pathogenic and likely pathogenic records used | 346,571 | printed by `scripts/differential_and_secondary_scan.py` (to standard error) |
| R11 | The 43-gene panel | BUB1B is the only panel gene in which the child carries a ClinVar pathogenic or likely pathogenic allele | `make expanded-differential` (written to `expanded-differential.txt` by `run_all.sh`) |
| R12 | ClinVar's BUB1B missense records (snapshot 2026-09-13) | 1,426 of uncertain significance, 35 benign or likely benign, 28 conflicting, 1 likely pathogenic (1,490 records; 96% uncertain). Pathogenic and likely pathogenic records by type: 37 nonsense, 35 frameshift, 10 splice site, 1 missense, 1 intronic | `bcftools query -r 15:40160984-40221137 -f '%ID\t%INFO/CLNSIG\t%INFO/MC\n' clinvar.vcf.gz` |
| R13 | gnomAD v4.1 counts for the two variants | stop-gain: exomes 115 of 1,461,846 alleles (AF 7.9e-5), genomes 5 of 152,174, never homozygous (120 alleles in total); missense: exomes 1 of 1,461,878, none in genomes. The exome frequencies also appear in `biallelic-ranking.txt` | gnomAD browser, variants `15-40209701-T-G` and `15-40220612-T-G` |
| R14 | Conservation of the two residues | In 115 vertebrate BubR1 sequences (UniProt orthologues aligned pairwise to human, not a curated multiple alignment): residue 737 is Leu in 110, Phe in 2, Ser in 2, Pro in 1; residue 1002 is Asn in 111, His in 1, gap in 3 | UniProt orthologue entries, taxonomy Vertebrata |
| R15 | Experimental structures with the same asparagine clamp | The related kinases Drosophila BubR1 (PDB 6JKK, 6JKM) and human BUB1 (PDB 6F7B, 4QPM, 5DMZ) | open the PDB entries and measure the contacts of the equivalent Asn |
| R16 | gnomAD variant co-occurrence service | Returns "not available for gnomAD v4"; on a control pair in an unrelated gene with nobody carrying both variants it returns `p_compound_heterozygous = 1`, which is the service having nothing to fit | a GraphQL `variant_cooccurrence` query on two public variants |
| R17 | Run time and hardware | About 55 minutes wall clock from the VCF to the predictions file on 4 vCPU, 16 GB RAM and 30 GB disk, most of it downloading references (run September 2026 (earlier package)) | `time ./run_all.sh` on such a machine |
| R18 | Reach of the gene's regulatory signals | Public GTEx maps list 8 independent expression-controlling (eQTL) signals for BUB1B between 0.78 Mb upstream and 0.84 Mb downstream of the transcription start site, and 153 splicing signals out to 1 Mb. The window we searched (61 kb before the gene start to 49 kb after its end) excludes 6 of the 8. A later GTEx v8 check (2026-10-25) found no Whole_Blood signal and only one testis signal among the eight (the strongest, p = 5e-17, inside the searched window, as is a signal in EBV-transformed lymphocytes, p = 1.2e-5); the six outside the window are in cerebellum (three), tibial artery, oesophagus mucosa and skin, with p-values from 4e-5 to 2e-10, all weaker than the testis signal. Tissue specificity was not tested | GTEx portal, gene BUB1B (the script is not part of this repository) |
| R19 | Aneuploidy screens | Allele-balance spread is higher on chromosomes 20, 21, 22 (and 9 and 17), attributed to repeats, low mappability and pericentromeric regions, with no two-peaked pattern; overdispersion of allele balance over 2,087,774 sites is 1.27 times the binomial expectation and tracks chromosome GC content (Spearman 0.72, p < 0.001); chromosome 20's excess is confined to 20 to 40 Mb. No copy-neutral loss of heterozygosity on 15q (median heterozygous fraction 0.068 against 0.067 genome-wide) (run 2026-09-14 and 2026-10-15) | `make overdispersion` and `make roh` |
| R20 | Protein language models and structure predictors on p.Asn1002Lys | 12 language models (ESM2 35M to 3B, ESM1b, ESM-1v five seeds, ESM-C) put the change between the 14th and 79th percentile of damage for this gene (median 41st); lysine is among the mildest changes at the site in every model. Four structure predictors (Boltz-2, OpenFold3, ESMFold v1, ESMFold2) return the unchanged fold and show only the lost hydrogen bonds; none separated labelled disease alleles or synthetic damaging controls from neutral ones, so none is used as evidence (run 2026-10-15 to 2026-10-16) | research code from the Track 2 work, not included here |
| R21 | The first mitochondrial pass with a fixed 2% cutoff | 1,044 variants: 10 homoplasmic and 1,034 low-level, two of which matched ClinVar pathogenic records. All were sequencing error or copies of mitochondrial sequence in the nucleus. The final run (`mitochondrial-screen.txt`) sets the floor from the data at 21.2% and counts 555 positions above 2% at its depth filter (run 2026-09-21) | `make mito-extract && make mito-call`, then lower the floor to 2% |

Six files in this folder were edited by hand after they were generated, to match the scripts
as they now stand, and no number in them changed: `insilico-panel.txt` had one paragraph about
AlphaMissense replaced with the text now printed by `scripts/insilico_panel.py`;
`phenotype-reverse-check.txt` has one rival-disease label, the marker beside the MVA1 row and
two sentences changed to the wording the script now prints; `completeness-audit.txt` and
`acmg-classification.txt` have one sentence each reworded to the scripts' current text;
`differential-exclusion.txt` has its conclusion line reworded (no gene is completed by
ClinVar-listed pathogenic variants alone); `biallelic-ranking.txt` has one comment reworded.

A further round of hand edits, made without the data, brought these files in line with
scripts whose text or logic changed. No result value changed in any of them (label numbers such as an OMIM id and a PMID
were added to `crosscheck.txt`), and each now says so in its own header:
`phenotype-reverse-check.txt` (the rows of diseases that are not MVA-related are withheld; the supplied HPO terms are no longer called the child's
features; the SLF2 and SMC5 rows are labelled Atelis syndrome 1 and 2, as OMIM and ClinVar
name them), `crosscheck.txt` (the same gene labels, the BUB1 label, the source column for
CENATAC, and the CEP192 reference, which is Guo 2024, *HGG Advances*, PMID 37981762),
`differential-exclusion.txt` (the BUB1 description, and the competing diagnoses other than the MVA-related genes, which are
counted and no longer named: the unnamed count went from 12 to 24, with each newly hidden gene counted in the category the
public view gives it), `insilico-panel.txt` (the control column
that shows the most damaging value was labelled "max", and for SIFT that value is the
minimum; and the Gln467His row, where the REVEL and VEST4 calls are BP4 Moderate and not
BP4 Supporting, because `Calibration.strength` walked the benign bounds weakest-first and
capped every benign call at Supporting) and `mitochondrial-screen.txt` (a sentence that said
no calls were heteroplasmic, which the counts above it did not support; the file now gives the 12 calls that remain
after the 5 seam artefacts, but not their split into homoplasmic and sub-homoplasmic, which the current script prints).
Section 8.3 of the report was changed to match the counts in that file. Regenerating these
files on the data machine replaces the hand edits.

## Known open items

* **VEST4.** The values used for the VEST4 column of `insilico-panel.txt`, `acmg_classify.py`
  and the report (0.433 for p.Asn1002Lys, 0.983 for p.Leu1012Pro, 0.231 for p.Gln467His) do
  not reproduce against dbNSFP 4.8a as served by MyVariant.info (downloaded 2025-01-17): it
  returned 0.885 and 0.841, 0.901 and 0.960, and 0.305 and 0.201 on 2026-10-25, for the
  MANE transcript and the second transcript. The other four predictors reproduce exactly.
  The cause was not determined (a release difference is one possibility; it was not tested).
  One thing to check:
  `load_vest` takes the position of the MANE transcript in the track's "all-ids" column and
  uses it to index the score list after dropping "." entries, which would misalign the index
  if any score were missing. `make insilico-refs` fetches the UCSC track without a version,
  so re-running it does not pin the source. The class of allele 2 does not depend on VEST4
  (3 points at most), but the text that says VEST4 leans benign does.
* **Optional private files.** `differential-exclusion.txt` counts the competing-diagnosis
  genes that are not named. That count (24 in the committed file) needs `data/private/extra_differential.json`
  (a JSON object of gene symbol to description), which holds every competing-diagnosis gene except the MVA-related
  ones and is not in this repository on purpose. The phenotype check can likewise track further diseases from
  `data/private/extra_watch_diseases.json` (a JSON object of OMIM id to label); the public view hides those rows.
  Without it the script checks only the genes listed in the script and prints a note in the
  output. The "8 terms" in `phenotype-reverse-check.txt` need both `data/private/case_terms.json`
  and `data/private/withheld_terms.json`. A gene symbol that is misspelt in `extra_differential.json`
  prints as "no damaging variants", because the script cannot tell a typo from a clean gene.

## To close this list

On the machine that holds the data, run `./run_all.sh`, run the targets for R5 to R9 and
R12, commit the new files and delete the matching rows above. Then compare the new
numbers with the ones the report quotes. If any differ, the report is wrong and must be
corrected before anything is released.

Also read the full output of `scripts/rank_biallelic_genes.py` (the committed file shows only the
top rows) and confirm that BUB1B is the only gene with a ClinVar pathogenic or likely pathogenic
variant in a complete recessive-compatible genotype. The report states this in the Summary and
in section 4.4.
