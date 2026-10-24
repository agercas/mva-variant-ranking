# Methods description — MVA variant prediction

---

## Model number (fill out this form once per model/approach submitted; up to 6)

Model 1 of 1.

## Please describe your model/approach in detail.

We looked for the cause in two independent ways and checked the result in several others. All code is in the public repository; each step below names its output file in results/.

1. Effect of every variant. All 4,740,790 variants that passed the sequencing lab's filters (PASS) were annotated with bcftools csq (Ensembl release 116 gene model), one chromosome's reference sequence at a time. This needs about 250 MB of extra disk and keeps all data on the machine. 15,965 variants are predicted damaging.
2. Rank genes. The ranking takes no gene, disease or phenotype as input; ClinVar pathogenic labels, which carry disease knowledge, add to the score. Keep variants that are rare (gnomAD v4.1 frequency 0.001 or lower; where we had no frequency, absence from dbSNP served as a proxy), well mapped (mapping quality 55 or more), covered by 15 to 100 reads, and, if heterozygous, carried by at least a quarter of the reads. Score each gene by whether it has a genotype that fits recessive inheritance (one homozygous damaging variant, or two heterozygous ones).
3. List every variant in BUB1B and in a 170 kb window around it, at any filter status. Of the 14 variants in the gene, two change the protein: a stop-gain and a missense. They form a compound heterozygous pair.
4. Test for a second variant that is hard to see: SpliceAI on the 15 variants called in the gene's region, exon depth from re-aligned raw reads, the upstream region.
5. Look for other explanations: all 346,571 ClinVar pathogenic and likely pathogenic records, named competing diagnoses, a 43-gene panel, chromosomes X and Y, the mitochondrial genome (rebuilt from the raw reads), and a whole-chromosome aneuploidy screen.
6. Check from the other side: from the disease to its known genes, and from the child's HPO terms to all diseases. MVA1 ranks 3rd of 12,867 diseases.
7. Verify: we re-derived both protein changes by hand from the reference with a hard-coded genetic code, and re-made both calls from the raw reads with a different aligner.
8. Classify both variants with the ACMG/AMP point system and show what would change each class.

## 🔴 (required) If commercially-available Generative AI was used, record the provider, the plan or tier, and the relevant setting in your methods description. A line is enough, for example, "Anthropic API, Claude Sonnet, commercial terms, no training on customer content."

Anthropic Claude (Claude Code), commercial terms, no training on customer content.

## Please state whether the submission file is the automated output of your computational approach (preferred), or if it has undergone downstream manual review and curation.

Both. The candidate variants are the automated output of the ranking, which takes no phenotype, disease name or gene as input. A script (scripts/build_track1_submission.py) writes the predictions file and checks it against the organizers' format and a copy of their scoring rules; the set of rows and their order were chosen by hand.

## Please describe any downstream manual review in detail.

We made three decisions by hand; the report explains them (sections 1, 5.3 and 6):

1. Rows 2 and 3 cover a second variant that this data cannot show (a balanced rearrangement, a regulatory variant). Their presence and order are a judgement about residual uncertainty.
2. The epcr column holds ranking scores that order the rows, best first. The format asks for a probability, which needs outcome data we lack.
3. We applied the ACMG classes by hand on the point scale. Every criterion that earns no points has a stated reason in results/acmg-classification.txt.

## Please state whether your approach only used publicly available data (preferred) or if proprietary data was also used. If proprietary data were used, where possible, please also include a submission based only on publicly available data.

Publicly available data only.

## Please describe any public data used in detail (e.g., reference panels, population databases, literature).

Reference data were downloaded in bulk and queried on the local machine. The gnomAD browser and Ensembl received the identifiers of the two reported variants and the secondary finding. NVIDIA NIM, Biohub ESM and ColabFold received a protein sequence and a small number of candidate sites, for the language-model and structure checks in report section 8.1 (R20).

- Ensembl release 116 (GRCh38 gene models and reference sequence); Ensembl REST for an independent transcript model and gene positions; HGNC for gene symbols.
- ClinVar GRCh38 VCF, snapshot 2026-09-13.
- gnomAD v4.1 exomes (and genomes for the two reported variants).
- HPO release 2026-09-02 (phenotype.hpoa, hp.obo).
- AlphaMissense (public per-protein tables, AlphaFold DB); REVEL v1.3; VEST4, PolyPhen-2 and SIFT through UCSC dbNSFP tracks; SpliceAI 1.3.1 run locally.
- UniProt O60566; AlphaFold model AF-O60566-F1 v6.
- ClinGen gene-disease validity and dosage curation for BUB1B; GTEx for the reach of the gene's regulatory signals.
- GRCh38 analysis set (NCBI) for re-alignment; rCRS for the mitochondrial genome.
- Literature: 13 papers, listed in the report (section 10).

Tools: bcftools 1.19 to 1.24, minimap2, bbmap 40.02, mappy 2.31, samtools, Biopython, Python 3.11.

## Please describe any proprietary data used in detail.

None.

## Is your approach able to output proposed pairs of compound heterozygous candidate variants, or only single candidate variants?

Yes. A compound heterozygous pair is one row of the predictions file, written with the chrom_2, pos_2, ref_2 and alt_2 fields. Pair rows sit at ranks 1 to 3, where the scorer awards the most points, and the file builder keeps each pair in one row.

Phase is undetermined (the variants lie 10,911 bases apart and the longest fragment at the locus is 800 bases), so ACMG criterion PM3 (variant in trans with a pathogenic allele) earns no points. We prefer the in-trans arrangement because the child is affected and no other cause was found.

## How did you handle secondary or incidental findings, if any? (See the finding_type column in your submission.)

One secondary finding is reported, in row 5 (finding_type secondary, epcr 0.1, last so that the score is unaffected): LZTR1 p.Tyr748Ter (chr22:20996720 C>G), a stop-gain, ClinVar 1409252 Pathogenic/Likely pathogenic, heterozygous. Heterozygous loss of LZTR1 function predisposes to schwannomatosis (dominant, incomplete penetrance), for which clinical follow-up exists. The child's condition is attributed to BUB1B (report section 7). The finding needs confirmation in an accredited laboratory and review by the treating team.

Four further heterozygous pathogenic or likely pathogenic variants in unrelated genes are not reported, because sharing health information beyond the finding that bears on the child's presentation is the family's decision.

## Please provide an estimate of run time and cost for your approach (if possible).

Cost: compute only. All data and tools are free and public, and the approach needs no commercial licence, cloud service or GPU.

Main path (VCF to predictions file), measured on 4 vCPU, 16 GB RAM and 30 GB disk: about 55 minutes wall clock, most of it downloading public references (about 25 minutes, 7.7 GB). The annotation step takes about 20 minutes and the ranking and checks about 8 minutes. The raw-read checks (re-alignment, mitochondrial genome) use the 85 GB of FASTQ files; the mitochondrial step streams them, stores nothing and takes about 30 minutes.

Scalability: the scan holds one chromosome's reference at a time, so memory and disk stay constant as the number of genomes grows. The ranking takes no gene, disease or phenotype as input (it uses ClinVar labels), so the same scripts run on any single-sample GRCh38 VCF.

## Provide a method abstract (up to 500 words). Please include strengths and limitations of your method as well as methodological detail.

We report a compound heterozygous genotype in BUB1B: chr15:40209701 T>G (p.Leu737Ter) and chr15:40220612 T>G (p.Asn1002Lys), consistent with mosaic variegated aneuploidy syndrome 1 (OMIM 257300, recessive).

Method. We annotated all 4,740,790 variants that passed the sequencing lab's filters (PASS) with bcftools csq, one chromosome's reference at a time. We ranked genes by whether they carry a rare (gnomAD frequency, with absence from dbSNP as a proxy where we had none), well-mapped genotype that fits recessive inheritance, with no phenotype, disease or gene as input. BUB1B ranks first of 175 genes, first of 27 with real gnomAD frequencies everywhere, and first at all 108 combinations of the filter thresholds. In the gene, two of 14 variants change the protein.

Strengths. The ranking takes no gene, disease or phenotype as input, so the method suits the next undiagnosed child; ClinVar labels carry the disease knowledge and identify the gene (see Limitations). Both changes were re-derived by hand and re-called from raw reads with a second aligner. The result was also reached from the disease side (1 of 9 known MVA genes carries a pathogenic ClinVar variant) and checked against the child's features (MVA1 ranks 3rd of 12,867 diseases; the two diseases above it have no damaging variant in the child). All of ClinVar, a 43-gene panel, chromosomes X and Y and the mitochondrial genome (rebuilt from raw reads with a noise floor measured from the data) give one explanation. Allele 1 is Likely pathogenic on PVS1 alone. Computational predictors flag allele 2 as damaging (AlphaMissense 0.92, likely pathogenic; PolyPhen-2 and SIFT agree), worth at most 2 ACMG points, or 3 with the rarity criterion PM2. Numbers that rest on runs whose output files are not yet committed are listed in results/RECORDED-CHECKS.md.

Limitations. The ranking needs ClinVar labels to identify the gene: with ClinVar removed, BUB1B ranks 41st of 175. Allele 2 is a variant of uncertain significance (1 point). A single calibrated predictor cannot lift it: the strongest band is Strong (+4), which gives 5 points against the 6 needed, and the five predictors we scored disagree. No laboratory test of allele 2 exists. Phase is undetermined with this library, so PM3 earns no points. A second hit that short reads cannot see (balanced rearrangement, repeat expansion, methylation change) or that lies beyond our 170 kb window remains possible. The aneuploidy screen sees gains in about 13% of cells, so the diagnosis rests on genotype. Mitochondrial calls reach down to about 21% of copies.

Allele 2 can reach Likely pathogenic with a test validated to the very strong level, or with two pieces together, for example phase at the strong level with a predictor or a laboratory test.
