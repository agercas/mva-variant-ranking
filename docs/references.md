# Tools, data and papers

Everything the analysis used or relied on, with a link. Versions are the ones we ran. The
report (section 10) has the same list in short form. Numbers such as R5 refer to
[`results/RECORDED-CHECKS.md`](../results/RECORDED-CHECKS.md).

## Software

| Tool | Version | Used for | Run by | Link |
|---|---|---|---|---|
| bcftools (`csq`) | 1.19 to 1.24 | predicts the effect of each variant on the protein | `coding_consequence_scan.sh` | https://github.com/samtools/bcftools |
| samtools | | reads, indexes and pileups | `targeted_realign.sh`, `mito_*.py` | https://github.com/samtools/samtools |
| minimap2 and mappy | mappy 2.31 | re-aligns the raw reads | `targeted_realign.sh` | https://github.com/lh3/minimap2 |
| BBMap | 40.02 | read processing for the re-alignment | `targeted_realign.sh` | https://sourceforge.net/projects/bbmap/ |
| SpliceAI | 1.3.1, run locally | predicts splice changes (R5) | command in R5 | https://github.com/Illumina/SpliceAI |
| Biopython | | reads the AlphaFold structure | `n1002_geometry.py` | https://biopython.org |
| Python, pytest, ruff | Python 3.11 | scripts, tests, lint | all | https://www.python.org |

## Data

| Data | Version | Used for | Link |
|---|---|---|---|
| The competition dataset (gated) | | the child's VCF and reads | https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data |
| Competition page and rules | | the task, format and scorer | https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026 |
| Ensembl | release 116, GRCh38 | gene models and reference sequence | https://www.ensembl.org |
| Ensembl REST | | a second gene-position and transcript source (`crosscheck.py`) | https://rest.ensembl.org |
| HGNC | | gene symbols and positions | https://www.genenames.org |
| ClinVar | snapshot 2026-09-13 | variant classifications | https://www.ncbi.nlm.nih.gov/clinvar/ |
| gnomAD | v4.1 | variant frequencies | https://gnomad.broadinstitute.org |
| HPO | release 2026-09-02 | disease and feature annotations | https://hpo.jax.org |
| OMIM | MVA1, 257300 | the disease entry | https://omim.org/entry/257300 |
| ClinGen, *BUB1B* | checked 2026-10-21 | gene-disease rating and dosage page | https://search.clinicalgenome.org/kb/genes/HGNC:1149 and https://search.clinicalgenome.org/kb/gene-dosage/HGNC:1149 |
| UniProt O60566 | | the protein, its domains | https://www.uniprot.org/uniprotkb/O60566 |
| AlphaFold model AF-O60566-F1 | v6 | structure of BubR1 | https://alphafold.ebi.ac.uk/entry/O60566 |
| AlphaMissense tables | public per-protein tables | missense scores (`alphamissense_bub1b.py`) | https://github.com/google-deepmind/alphamissense |
| REVEL | v1.3 | missense score | https://sites.google.com/site/revelgenomics/ |
| VEST4, PolyPhen-2, SIFT | through UCSC dbNSFP tracks | missense scores | https://genome.ucsc.edu |
| GTEx | | reach of the gene's regulatory signals | https://gtexportal.org |
| rCRS and the GRCh38 analysis set | | mitochondrial reference; re-alignment reference | https://www.ncbi.nlm.nih.gov/nuccore/NC_012920 |

Hosted prediction and structure services were used once, for the checks in R20. They are
described there and are not part of the public code.

## Papers

| # | Paper | Access | Link |
|---|---|---|---|
| 1 | Hanks 2004, *Nat Genet* 36:1159. *BUB1B* mutations, constitutional aneuploidy and cancer | Paywalled | https://pubmed.ncbi.nlm.nih.gov/15475955/ |
| 2 | Suijkerbuijk 2010, *Cancer Res* 70:4891. Causes of BUBR1 dysfunction in MVA | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/20516114/ |
| 3 | Ochiai 2014, *PNAS* 111:1461. An intergenic variant upstream of *BUB1B* | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/24344301/ |
| 4 | Richards 2015, *Genet Med* 17:405. ACMG/AMP variant standards | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/25741868/ |
| 5 | Tavtigian 2020, *Hum Mutat* 41:1734. The ACMG/AMP point system | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/32720330/ |
| 6 | Abou Tayoun 2018, *Hum Mutat* 39:1517. The PVS1 criterion | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/30192042/ |
| 7 | Pejaver 2022, *Am J Hum Genet* 109:2163. Calibration of missense predictors | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/36413997/ |
| 8 | Bergquist 2025, *Genet Med* 27:101402. Calibration of further predictors, including AlphaMissense | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/40084623/ |
| 9 | Brnich 2019, *Genome Med* 12:3. The PS3/BS3 criterion | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/31892348/ |
| 10 | Piotrowski 2014, *Nat Genet* 46:182. *LZTR1* and schwannomatosis | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/24362817/ |
| 11 | Yamamoto 2015, *J Med Genet* 52:413. *LZTR1* and Noonan syndrome | Paywalled | https://pubmed.ncbi.nlm.nih.gov/25795793/ |
| 12 | Johnston 2018, *Genet Med* 20:1175. Recessive Noonan syndrome and *LZTR1* | Free (PMC) | https://pubmed.ncbi.nlm.nih.gov/29469822/ |
| 13 | Bohers 2008, *Hum Genet* 124:473. Dosage dependence of BUBR1 in HeLa cells | Paywalled | https://pubmed.ncbi.nlm.nih.gov/18932004/ |

Access: *Free (PMC)* is open to anyone on PubMed Central; *Paywalled* needs a subscription.

Papers describing the tools:

| Tool | Paper | Link |
|---|---|---|
| bcftools csq | Danecek & McCarthy 2017, *Bioinformatics* 33:2037 | https://pubmed.ncbi.nlm.nih.gov/28205675/ |
| minimap2 | Li 2018, *Bioinformatics* 34:3094 | https://pubmed.ncbi.nlm.nih.gov/29750242/ |
| SpliceAI | Jaganathan 2019, *Cell* 176:535 | https://pubmed.ncbi.nlm.nih.gov/30661751/ |
| gnomAD v4 | Chen 2024, *Nature* 625:92 | https://pubmed.ncbi.nlm.nih.gov/38057664/ |
| HPO | Köhler 2021, *Nucleic Acids Res* 49:D1207 | https://pubmed.ncbi.nlm.nih.gov/33264411/ |
| AlphaMissense | Cheng 2023, *Science* 381:eadg7492 | https://pubmed.ncbi.nlm.nih.gov/37733863/ |
| AlphaFold | Jumper 2021, *Nature* 596:583 | https://pubmed.ncbi.nlm.nih.gov/34265844/ |
| REVEL | Ioannidis 2016, *Am J Hum Genet* 99:877 | https://pubmed.ncbi.nlm.nih.gov/27666373/ |
| VEST4 | Carter 2013, *BMC Genomics* 14 (Suppl 3):S3 (the VEST method) | https://pubmed.ncbi.nlm.nih.gov/23819870/ |
| PolyPhen-2 | Adzhubei 2010, *Nat Methods* 7:248 | https://pubmed.ncbi.nlm.nih.gov/20354512/ |
| SIFT | Ng & Henikoff 2003, *Nucleic Acids Res* 31:3812 | https://pubmed.ncbi.nlm.nih.gov/12824425/ |
