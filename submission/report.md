# Report: a compound heterozygous *BUB1B* genotype in PROBAND01

---

## Summary

**Answer.** The child carries two different changes in the gene *BUB1B*. The gene encodes a
protein, BubR1, that holds a dividing cell until every chromosome is correctly attached. With
too little BubR1, cells divide with the wrong number of chromosomes. The matching condition is
mosaic variegated aneuploidy syndrome 1 (MVA1; OMIM 257300, an entry in the Online Mendelian
Inheritance in Man catalogue). It is recessive: both copies of the gene must be affected. A
different damaging change on each copy is called compound heterozygous.

The table gives the ACMG class of each change. ACMG is the standard grading of variants into
Pathogenic, Likely pathogenic, VUS (uncertain significance), Likely benign and Benign.

| | Variant (GRCh38) | Effect | ACMG class |
|---|---|---|---|
| Allele 1 | `chr15:40209701 T>G`, `p.Leu737Ter` | stop signal in the middle of the gene | Likely pathogenic |
| Allele 2 | `chr15:40220612 T>G`, `p.Asn1002Lys` | one amino acid swapped, in the last exon | Uncertain significance (VUS) |

(Transcript `NM_001211.6`: `c.2210T>G` and `c.3006T>G`.)

**Confidence.** We are confident in the gene and in allele 1. Allele 2 stays a variant of
uncertain significance. The `epcr` column (the organizers' Estimated Probability of Causal
Relationship) holds ranking scores, and the competition's scorer reads only the order of the
rows (sections 1 and 6.1).

**The evidence.**

* ClinVar (a public database of variant classifications) lists allele 1 as pathogenic or
  likely pathogenic for MVA1.
* We ranked the genes of the whole genome by whether they carry a genotype that fits
  recessive inheritance. *BUB1B* came first of 175 genes. The ranking uses ClinVar labels;
  without the label on allele 1, *BUB1B* is 41st (section 4.1).
* Two of the 14 variants in *BUB1B* change the protein, and they form one pair (section 4.3).
  *BUB1B* is the only gene where a ClinVar-pathogenic variant is part of a complete recessive
  genotype (section 4.4). Starting from the disease leads to the same gene (section 4.6).
* A stop on one copy plus an amino-acid swap on the other is the pattern published for
  MVA1 patients (section 5.2).
* Computational predictors flag allele 2 as damaging: AlphaMissense scores it 0.9229 (likely
  pathogenic), and PolyPhen-2 and SIFT agree. REVEL gives no evidence and VEST4 leans benign.
  For these scores ClinGen's calibration gives one predictor at most 2 points. With 1 point for
  rarity that makes 3, and Likely pathogenic needs 6, so the formal grade stays VUS (sections 5.3
  and 8.1).

**Open points.**

* No laboratory has tested allele 2, and prediction programs disagree about it, so it stays
  "uncertain" (section 5.3).
* Phase (which variants share a copy) is undetermined: the sequenced DNA fragments reach 800
  bases and the variants lie 10,911 bases apart (section 6.2).
* A second variant that short reads miss, or that lies far from the gene, remains possible.
  Rows 2 and 3 of the predictions cover this.

Numbers marked R1 to R21 come from recorded runs whose output files are not yet in the
repository (section 3).

**Contents.**

1. The predictions and how they are scored
2. The condition, in plain words
3. Method
4. Results
5. Interpretation
6. Confidence, and what could be wrong
7. Secondary finding
8. Additional findings
9. Reproducibility, scalability, privacy and use of AI
10. Tools, data and references

---

## 1. The predictions and how they are scored

The file has five rows, best first. The `epcr` column, which the format describes as a
probability, holds a ranking score: it orders the rows from most to least likely, and the scorer
reads only that order.

| Row | `epcr` (ranking score) | Content | Role of the row |
|---|---|---|---|
| 1 | 0.9 | Both variants, in one row | Our answer |
| 2 | 0.6 | Allele 1 alone | Covers a second variant that this data cannot show |
| 3 | 0.4 | Allele 1 plus `chr15:40216470 A>G`, a change far inside an intron | The only other novel (absent from dbSNP) heterozygous variant in the gene. A weak candidate: SpliceAI (a program that predicts splicing changes) predicts no effect |
| 4 | 0.2 | Allele 2 alone | Covers allele 2 as the sole answer |
| 5 | 0.1 | A secondary finding in *LZTR1* (section 7) | The organizers invited secondary findings. It goes last, and F-max takes the best cut-off, so the row leaves the score unchanged |

The scorer (reimplemented in `src/mva/scoring.py` from the organizers' published evaluation code):

| Rank of the row that matches the true answer fully | 1 | 2 or 3 | 4 or 5 | 6 to 10 |
|---|---|---|---|---|
| Rank points | 100 | 50 | 25 | 10 |

* It looks for a row that matches the true answer exactly, and stops at the first one. If
  there is none, one correct variant of a pair earns half the points of its row.
* A pair counts only if both variants are in one row. Two rows with one variant each give
  two half matches.
* A pair at rank 4 or lower scores 25 or 10, below the 50 that half credit gives, so rows that hold
  a pair sit at ranks 1 to 3.
* F-max is the best balance of precision and recall over all cut-offs on `epcr`, counted per
  variant. It ignores `finding_type`, so the secondary finding takes a rank.

Scores of the file under three possible true answers, computed with a copy of these rules:

| If the true answer is | Rank points | F-max |
|---|---|---|
| our pair | 100.0 | 1.000 |
| allele 1 plus the intronic variant | 50.0 | 0.800 |
| allele 1 alone | 50.0 | 0.667 |

Row 1 holds both variants, so it earns half credit when the true pair shares one variant with it.
Row 3 raises F-max in the second case (0.800 instead of 0.500). Row 2 earns the 50 points in the
third case; without it that case scores 0. The cost is small: if allele 2 alone were the truth, row 2 pushes it from
rank 3 (50 points) to rank 4 (25).

---

## 2. The condition, in plain words

Before a cell divides it copies every chromosome. A checkpoint then waits until each copy is
attached correctly to the machinery that pulls the two sets apart. BubR1 is part of this
checkpoint. With too little BubR1 the checkpoint is weak, cells divide early, and the daughter
cells end up with too many or too few chromosomes (aneuploidy). In mosaic variegated aneuploidy
(MVA) a share of the body's cells is affected, each with different errors. Children with MVA1
have growth restriction and a high risk of childhood cancer. The competition describes fewer
than 50 known patients worldwide.

Each person has two copies of every gene. A recessive condition needs both copies damaged;
*compound heterozygous* means the two copies carry two different damaging changes.

| Word | Meaning here |
|---|---|
| variant, allele | a difference from the reference sequence; one copy of it |
| stop-gain (nonsense) | a change that creates a stop signal in the middle of a gene |
| missense | a change that swaps one amino acid in the protein for another |
| NMD | nonsense-mediated decay: the cell destroys messages that contain an early stop signal. A stop more than 50 bases before the last exon-exon junction (where the final two exons join) is expected to trigger it |
| carrier | a healthy person with one damaged copy of a recessive gene |
| in *trans* / in *cis* | on different copies of the gene / on the same copy |
| phase | which variants lie on the same copy |
| VUS | variant of uncertain significance |
| ACMG/AMP points | each piece of evidence earns points (1 supporting, 2 moderate, 4 strong, 8 very strong; Tavtigian 2020). 10 or more is Pathogenic, 6 to 9 Likely pathogenic, 0 to 5 VUS |
| ClinVar, gnomAD | public databases of variant classifications, and of variant frequencies in about 800,000 people |
| dbSNP | a public catalogue of known variants; "not in dbSNP" serves as a rough sign of rarity |
| VCF, FASTQ | the sequencing lab's table of the child's variants; the raw sequencing reads |
| indel | a small insertion or deletion |
| ClinGen | the Clinical Genome Resource, which rates gene–disease links and issues guidance for classifying variants |
| GTEx | public maps of how variants change gene expression in body tissues |
| shRNA | a short RNA that lowers the amount of one protein in cultured cells |
| HPO | the Human Phenotype Ontology, a standard list of clinical features |

---

## 3. Method

We worked from two starting points: the genome and the disease. The genome-wide ranking takes
no gene, disease or phenotype as input; its score uses ClinVar labels, which carry disease
knowledge (4.1). The disease-first cross-check starts from MVA. Disclosure: early in the work we
saw file names on the public leaderboard that name *BUB1B*. Every row of our file and its
position come from our own analysis of the genome.

1. **Start.** The 4,740,790 variants that passed the sequencing lab's filters (R1).
2. **Effect on the protein.** We predicted the effect of each variant with `bcftools csq` and
   the Ensembl gene model, one chromosome at a time to keep disk use low (about 250 MB extra).
   15,965 variants are predicted damaging.
3. **Filters.** We kept variants that are rare (gnomAD frequency 0.001 or lower; where we had
   no frequency, "not in dbSNP" served as a proxy, see 4.1), well mapped (mapping quality 55 or
   more), covered by 15 to 100 reads and, if heterozygous, carried by at least a quarter of the
   reads. These filters remove false calls, where reads from repeated sequence land in the wrong
   place and crowd the top of the list.
4. **Ranking.** Genes score higher if they carry a homozygous damaging variant or two
   heterozygous ones, and higher again if ClinVar calls one of the variants pathogenic.
5. **Variants in the gene.** We listed every variant in the top gene and in a 170 kb window
   around it, including variants that failed the lab's filters.
6. **A second variant that is hard to see.** We tested for a splice change, an exon deletion
   and a variant just upstream of the gene (R5, R6, R9).
7. **Other explanations.** We checked all ClinVar pathogenic records, named competing
   diagnoses, a 43-gene panel, chromosomes X and Y, the mitochondrial genome and a
   whole-chromosome aneuploidy screen.
8. **The other direction.** We worked from the disease to its known genes, and from the
   child's features to all diseases.
9. **Verification.** We re-derived the effect of each variant by hand and re-called both
   variants from the raw reads with a different aligner (R7).
10. **Classification.** We graded both variants with the ACMG rules and asked what would
    change each grade.

`MANIFEST.md` in the repository lists the script for each step and the file that holds its
output. R1 to R21 are entries in `results/RECORDED-CHECKS.md`: runs we made on the machine that
holds the data and wrote down at the time, with output files not yet in `results/`. A reader can
check those numbers once the files are added.

### 3.1 The path from variants to answer, step by step

The numbers at each stage (from `biallelic-ranking.txt`, `crosscheck.txt` and R1 to R9):

```
 4,740,790  variants that passed the lab's filters
     |      step 2: predict the effect on the protein (bcftools csq)
    15,965  variants predicted damaging
     |      step 3: drop common (13,270 by dbSNP proxy, 386 by gnomAD), badly mapped (1,179),
     |              lopsided allele balance (273), wrong depth (148), unmatched (27)
       682  rare, well-mapped damaging variants
     |      step 4: keep genes whose genotype could cause a recessive disease
       175  genes      <- BUB1B ranks 1st (41st if ClinVar labels are removed)
     |      redo with real gnomAD frequencies on every chromosome
        27  genes      <- BUB1B ranks 1st; also 1st at all 108 filter settings (14 to 64 genes)
     |      steps 5 to 7: look in the gene, look for a second variant, look for other explanations
         1  gene, 2 variants (stop-gain + missense), one explanation
```

The work went through four stages, in this order. The first three each changed an earlier belief.

| Stage | What we did | What we found | What changed |
|---|---|---|---|
| Genome-wide ranking | Steps 1 to 4 with the initial filter thresholds, using "no dbSNP ID" as a proxy for rare | *BUB1B* first of 209 genes, with a stop-gain and a missense variant. The ranking in `biallelic-ranking.txt`, which adds the annotation fix and the allele-balance filter described below but keeps the proxy, gives 175 genes | The proxy agreed with real gnomAD frequencies in only 35.7% of cases on chromosome 15 |
| Corrected ranking | Real gnomAD frequencies, a fix for an annotation bug (9.1% of records were never protein-changing) and an allele-balance filter | The candidate pool shrank from 209 to 48, 40 and then 27 genes; *BUB1B* stayed first at all 108 threshold settings | The ranking depends on ClinVar (4.1) |
| Search of the gene | Variants in a 170 kb window, splice scores, exon depth, upstream regulatory region, re-calling from raw reads | Exon depth, splice scores and the upstream region all looked normal. Both calls held up (allele fractions 0.543 and 0.517) | Phase is undetermined with this library (6.2) |
| Search for other explanations | Other pathogenic variants genome-wide, named competing diagnoses, a 43-gene panel, X, Y, mitochondria, aneuploidy; working from the disease and from the child's features | One explanation. The features fit MVA1 and fit other diseases too | None; the answer stood |

`MANIFEST.md` lists the script and output file for each step.

---

## 4. Results

### 4.1 The top-ranked gene, and its dependence on ClinVar

| Run | *BUB1B* |
|---|---|
| Ranking with "not in dbSNP" as a proxy for rare outside chromosome 15 | first of **175** genes |
| Re-run with a real gnomAD frequency for every variant (R2) | first of **27** genes |
| That re-run at 108 combinations of the filter thresholds (R3) | first in **all 108** |
| The same ranking with ClinVar labels removed | **41st** of 175 |

The genotype pattern narrows 4.7 million variants to 175 candidate genes. The ranking needs the
ClinVar label on allele 1, which names MVA1, to put *BUB1B* first.

### 4.2 The two variants

| | Allele 1 | Allele 2 |
|---|---|---|
| Position (GRCh38) | `chr15:40209701 T>G` | `chr15:40220612 T>G` |
| Protein change | p.Leu737Ter, codon TTA to TGA | p.Asn1002Lys, codon AAT to AAG |
| Location | exon 17 of 23 | exon 23 of 23 (last exon) |
| Distance to the last exon-exon junction | 747 bases | inside the last exon |
| gnomAD v4.1 | 115 of 1,461,846 exome alleles (AF 7.9e-5) and 5 of 152,174 genome alleles, never homozygous (R13) | 1 of 1,461,878 exome alleles, none in genomes (R13) |
| ClinVar | 533901, Pathogenic/Likely pathogenic for MVA1 | not recorded as this nucleotide change. ClinVar 4600147 (`c.3006T>A`) gives the same protein change and is a VUS |
| Reads that carry the variant | 25 of 46 | 13 of 28 |
| Alternate-allele fraction when both are re-called from raw reads with a second aligner (R7) | 0.543 | 0.517 |

Allele 1 is seen 120 times in about 1.6 million gnomAD alleles, a carrier frequency of about 1
in 6,700.

We re-derived both protein changes by hand from the reference sequence and the standard codon
table. A control variant inside an intron was correctly reported as outside the coding
sequence (`consequence-verification.txt`).

### 4.3 Anything else in the gene

The search window runs from 61 kb before the gene to 49 kb after its end, in short-read data.
Short reads detect substitutions, small insertions and deletions, and exon-scale changes in
depth. Balanced rearrangements, deletions wholly inside an intron, repeat expansions and
methylation changes need other methods. The gene's regulatory reach extends beyond the window:
public GTEx maps of variants that control expression list signals for *BUB1B* up to 0.8 Mb
away, and six of the eight independent signals lie outside the window (R18). The child's
genotypes there are unexamined. "Two protein-changing variants" therefore describes this window
and what short reads see.

The gene has 14 variants, whether or not they passed the filters. Two change the protein, the
two above, and the other twelve lie outside the coding sequence. The 170 kb window holds 185
variants, and none of the others changes a BUB1B protein (`completeness-audit.txt`).

| Candidate second hit | Test | Result |
|---|---|---|
| a change that disrupts splicing | SpliceAI, run locally on the 15 variants called in the gene's region (R5) | all 0.03 or lower; a synthetic splice-site control scores 0.98 |
| a deletion of one or more exons | 22,234 read pairs from the raw reads, re-aligned (R6) | all 23 coding exons at 0.77 to 1.44 times the median depth; none near half |
| a regulatory change just upstream (a variant 44 kb upstream of *BUB1B* is the proposed cause of MVA in Japanese families, and edited cells reproduce the cellular features, Ochiai 2014) | the 61 kb before the gene (R9) | 88 variants; the 7 without a dbSNP ID are indels in repeats with the pattern of sequencing errors |

A deletion of a whole exon on one copy shows as a depth near 0.5. The depth test covers whole
exons; smaller deletions inside an exon need other methods.

### 4.4 Other explanations

* **All of ClinVar.** We compared every variant the child carries with all 346,571
  pathogenic or likely pathogenic ClinVar records (R10), without a gene list. The child
  carries six such variants with good read support, all heterozygous, and no gene holds two of
  them (`differential-exclusion.txt`). *BUB1B* is the only gene where one of them is part of a
  complete recessive genotype, with allele 2, which ClinVar does not list as pathogenic
  (`biallelic-ranking.txt`, section 5.3). Section 7 describes the other five. The sweep covers
  variants ClinVar already lists.
* **Competing diagnoses.** The genes for the two diseases that rank above MVA1 on phenotype,
  and the genes for the other conditions we considered, hold no damaging variant or only variants
  that ClinVar calls benign or likely benign. *LZTR1* is the exception (section 7). One of the
  other genes carries a variant that ClinVar has not classified. The other known MVA genes
  (*CEP57*, *TRIP13*, *BUB1*, *CENATAC*, *SLF2*, *SMC5*) hold no ClinVar pathogenic variant
  (`crosscheck.txt`). Imprinting changes need methylation testing.
* **A 43-gene panel** of MVA, chromosome-breakage, Fanconi pathway and cancer-predisposition
  genes: *BUB1B* is the only gene in which the child carries a ClinVar pathogenic or likely
  pathogenic variant (R11). The panel covers variants ClinVar lists.
* **Sex chromosomes.** The child is male, so single-copy X variants appear as homozygous and
  the ranking covers X-linked recessive inheritance. The only two damaging variants on Y failed
  the call filters and sit in a gene unrelated to this condition (`completeness-audit.txt`).
* **Mitochondrial DNA.** The call set lacks the mitochondrial genome, so we rebuilt it from 273
  million raw reads (section 8.3). Calls reach down to variants in about 21% of mitochondrial
  copies, and none of the pathogenic or likely pathogenic variants is present at that level.
* **Whole-chromosome aneuploidy.** We screened allele balance and depth
  (`aneuploidy-screen.txt`). The screen detects a gain in about 13% of cells or more (section
  8.2), a level at which a mosaic condition stays undetected, so its result is uninformative
  here. Three chromosomes give estimates above 13% (21: 0.150, 22: 0.177, Y: 0.416). A male has
  one copy of Y, so heterozygous calls there are mapping artefacts. Chromosomes 21 and 22 are
  rich in repeats, which distort allele balance, and the two-peaked pattern of a real mosaic
  gain is absent from both (R19). We attribute these values to artefacts; that is an inference
  that direct tests have yet to confirm.

### 4.5 The child's clinical features and MVA1

The dataset comes with 8 HPO terms. The numbers below use all 8 as supplied.

* We ranked all 12,867 diseases in HPO by how well their annotated features match the terms,
  giving more weight to rare features. MVA1 ranks **3rd**. Two diseases rank
  above it, and the child carries no damaging variant in either gene.
* We graded each link between a feature and MVA1 by how specific the connecting HPO term is,
  because every term shares the ontology root with every disease. Two features are annotated on
  MVA1 directly. One more connects through a specific term, and two connect through general
  terms. MVA1 therefore accounts for 5 of the 8 features,
  3 of them informatively. The other three are very general features, shared with most diseases
  (`phenotype-reverse-check.txt`).

The child's features are compatible with MVA1 and support it weakly. The three very general
features need no separate explanation.

### 4.6 Starting from the disease

We also started from the disease, with different sources from the main analysis: a published
gene list for MVA (nine genes, from OMIM and public summaries), positions from HGNC and
Ensembl, and the child's variants read straight from the VCF. The child carries a ClinVar
pathogenic variant in one of the nine genes: *BUB1B* (`crosscheck.txt`). This route finds
the same gene and the same two coding variants. It starts from the disease, which leaves the
diagnosis itself untested.

---

## 5. Interpretation

### 5.1 Effects of the two changes

* **Allele 1** puts a stop signal at codon 737, 747 bases before the last exon-exon junction.
  Any stop more than 50 bases before that junction normally makes the cell destroy the
  message, so this copy should make little or no BubR1. Most pathogenic *BUB1B* variants in
  ClinVar are of this kind: 37 nonsense, 35 frameshift and 10 splice-site, and one missense
  (R12).
* **Allele 2** is in the last exon, so its message is not expected to be destroyed. The cell
  makes a full-length protein with lysine in place of asparagine at position 1002. In the
  AlphaFold model of BubR1 this asparagine sits mostly inside the protein (relative solvent
  accessibility 0.14) and makes three hydrogen bonds to backbone atoms, two to tryptophan
  978 and one to valine 998 (`n1002-geometry.txt`). The position holds asparagine in 111 of
  115 vertebrate sequences (R14). The same arrangement appears in solved structures of two
  related kinases (R15). Lysine has no side-chain oxygen, so it cannot accept the hydrogen
  bond that asparagine takes from tryptophan 978. This is a model-based prediction; the protein
  made from this allele is unmeasured. Residue 1002 lies in BubR1's kinase domain
  (UniProt residues 766-1050). Most authors regard it as a pseudokinase that works without enzyme
  activity, although this is not settled. The analysis therefore rests on hydrogen bonding and
  conservation.

If the two variants lie on different copies, nearly all of the BubR1 in the child's cells
comes from allele 2. For the condition to appear, that protein must be low in amount or poor
in function.

### 5.2 A stop plus a missense, and MVA1

The first patients linked to *BUB1B* came from five families in which truncating and missense
variants were found (Hanks 2004). In patient cell lines, Suijkerbuijk and colleagues (2010)
found three things.

* In patients with biallelic mutations, a missense variant near the kinase domain was paired
  with a truncating variant on the other copy.
* BubR1 protein was two to six times lower than normal, because the missense protein is
  turned over faster.
* Two of the low-abundance missense proteins (p.Ile909Thr and p.Leu1012Pro) worked normally
  when forced to high levels. Their defect was a low amount of protein, with function intact.

MVA1 patients are expected to have at least one copy that still makes some protein. A missense
change that lowers the protein level fits this; a second stop would be a surprise.

The amount of BubR1 in the child is unknown. Allele 1 gives about none, and allele 2 gives an
unmeasured fraction of what one copy normally makes.

Bohers and colleagues (2008) lowered BubR1 with shRNA in five stable HeLa clones that kept 8.5,
10, 14, 58 and 77% of normal BubR1 (Figure 1b).

* Premature chromatid separation rose as BubR1 fell, from 5% of control cells to 12% at 77% and
  39% at 8.5% (Figure 3b).
* The abstract states that aneuploidy was clearly detected below 50% residual BubR1, and
  ClinGen's dosage curation of *BUB1B* (last evaluated 2022) repeats that figure.
* The paper defines aneuploidy as more than 4N DNA after 48 to 72 hours in nocodazole (a drug
  that arrests mitosis), a measure of cells that escape the arrest and copy their DNA again
  (Figure 4b). The 58% clone matched controls at 48 hours (6 to 8% of cells) and began to rise at
  72 hours, and the 14% clone reached 18% at 48 hours. No clone lies near 50%, so the data place
  the change between 58% and 14%.
* The series comes from one cancer-derived cell line, selected clones and no patient cells. The
  BubR1 level at which a person develops disease is unknown.

Carriers of one damaged copy, such as the parents of patients, have about half the normal amount
and are generally described as unaffected. ClinGen's Hereditary Cancer expert panel rates the
*BUB1B*–MVA1 link Definitive, with recessive inheritance (2019).

The diagnosis depends on allele 2 being a hypomorph (an allele that works partly), a claim this
data cannot test.

### 5.3 The ACMG grades

Grades follow the ACMG/AMP standard (Richards 2015) on the point scale (Tavtigian 2020).

| | Allele 1 | Allele 2 |
|---|---|---|
| Counted | PVS1 (very strong, +8; Abou Tayoun 2018): a stop expected to switch the gene off, in a gene where patients are known to have lost its function | PM2 (supporting, +1): seen once in 1.46 million alleles |
| Considered, no points | PM2 (seen 120 times, so common enough to fail the criterion), PM3, PS3, PP5 | PP3, PM1, PM3, PS3, PP4 |
| Points | 8 | 1 |
| Grade | **Likely pathogenic** | **VUS** |

The codes: PM2 is absent or extremely rare in population databases. PM3 is a variant found in
*trans* with a pathogenic variant, which needs phase. PS3 is a laboratory test that shows
damage. PP5 is a reputable source's classification, which ClinGen no longer recommends as
evidence. PP3 is computer predictors agreeing that the variant is damaging. PM1 is a variant in
an important functional site or hot spot. PP4 is a phenotype specific to one gene; several genes
cause MVA.

Allele 1 is Likely pathogenic on PVS1 alone; adding PM2 would make 9 points and the same grade.
PM1 earns no points for allele 2: the residue is conserved and probably structural, and no
laboratory result marks it as critical.

**Computational predictors flag allele 2 as damaging.** AlphaMissense (a deep-learning model
built on AlphaFold) scores p.Asn1002Lys 0.9229 and labels it likely pathogenic. PolyPhen-2
(0.997) and SIFT (0.000) agree. REVEL (0.472) gives no evidence and VEST4 (0.433) leans benign.
The structure points the same way (5.1). Section 8.1 examines the AlphaMissense score.

ClinGen turns each score into ACMG points through calibrated bands (Pejaver 2022 for REVEL,
VEST4, PolyPhen-2 HumVar and SIFT; Bergquist 2025 for AlphaMissense) and asks for one predictor
chosen before the variant is seen. We report all five. The strongest bands for this variant give moderate evidence
(+2), so the best case reaches 3 points against the 6 that Likely pathogenic needs, and the formal
grade stays VUS:

| Predictor | Score for p.Asn1002Lys | Evidence under the ClinGen bands | Total points for allele 2 |
|---|---|---|---|
| REVEL | 0.472 | none | 1 |
| VEST4 | 0.433 | against pathogenic, supporting (−1) | 0 |
| PolyPhen-2 | 0.997 | for pathogenic, supporting (+1) | 2 |
| SIFT | 0.000 | for pathogenic, moderate (+2) | 3 |
| AlphaMissense | 0.9229 | for pathogenic, moderate (+2) | 3 |

Among the five, the highest total is 3 points and Likely pathogenic starts at 6
(`acmg-classification.txt`). The strongest band any calibrated predictor reaches is Strong (+4),
which gives 5 points, so a predictor alone cannot lift the allele. ClinGen has no
*BUB1B*-specific variant curation panel (checked 2026-10-21), so the general bands apply. For a well-studied MVA1
missense variant (p.Leu1012Pro, ten residues from ours) REVEL also gives no criterion
(`insilico-panel.txt`).

**Evidence that would classify allele 2.** It has 1 point and Likely pathogenic starts at 6:

| Added evidence | Total | Grade |
|---|---|---|
| a predictor chosen in advance, supporting / moderate / strong | 2 / 3 / 5 | VUS |
| phase in *trans* (PM3 moderate, +2) | 3 | VUS |
| phase in *trans* plus a second unrelated patient (PM3 strong, +4) | 5 | VUS |
| a test showing less BubR1 protein (PS3 strong, +4) | 5 | VUS |
| a test validated to the very strong level (PS3 very strong, +8; Brnich 2019) | 9 | Likely pathogenic |
| phase (strong) plus a predictor (supporting), or phase (strong) plus a test (strong) | 6 or 9 | Likely pathogenic |

(PM3 uses the 2015 default weights.) Phase is the most informative single experiment and needs
parental samples or long reads. Reclassifying allele 2 therefore requires new data; reanalysis of
the supplied files reaches at most 3 points.

---

## 6. Confidence, and what could be wrong

### 6.1 Confidence

A calibrated probability needs outcome data, such as how often calls like ours prove right,
and that data is unavailable. The table gives our confidence in words, with reasons.

| Statement | Confidence | Basis |
|---|---|---|
| The causal gene is *BUB1B* | high | ClinVar label on allele 1, the ClinGen rating Definitive, no competing diagnosis found, features compatible. The ranking alone would not give this (section 4.1) |
| Allele 1 is the stop-gain | high | ClinVar pathogenic, Likely pathogenic on PVS1, re-derived by hand and re-called from raw reads |
| Allele 2 is the missense | moderate | it is the only other protein change in the gene, the published pattern is a stop plus a missense, and AlphaMissense, PolyPhen-2 and SIFT flag it as damaging. It is formally a VUS. A second variant that short reads miss, a harmless missense, or both variants on one copy would each change the answer |

### 6.2 Phase

Phase is undetermined. Read-backed phasing joins variants close enough to appear on one DNA
fragment. The two variants are 10,911 bases apart and the longest fragment seen at this site
is 800 bases (median 437, R8). Every gap between neighbouring heterozygous sites at the site
exceeds 800 bases, so no chain of nearby variants links them (R8). Answering this takes a
different library or parental samples.

Two rare variants in one gene, with no other information, are about equally likely to lie on
the same copy or on different copies. The child's condition moves this towards different
copies, because the organizers state that the child has MVA and MVA1 is recessive. The table
shows the reasoning.

| Arrangement | What it means | Our view |
|---|---|---|
| In *trans* (one on each copy) | the expected compound heterozygote | most likely |
| Both on one copy, the other copy normal | the child would have about half the normal BubR1, like a carrier, and the condition would be unexplained | less likely; cannot be excluded without the parents |
| Both on one copy, plus a third change on the other copy that this data does not show | the condition is explained, but the reported pair is wrong | possible; rows 2 and 3 cover it |

This reasoning earns no ACMG points for PM3, which needs independent evidence of the *trans*
arrangement: arguing from "the child is affected, so the pair is causal" would be circular,
because it uses the conclusion as its own proof.

We also tried the gnomAD co-occurrence service, the one public service that estimates phase
without parents. It works only on v2 exomes, and for a control pair in an unrelated gene that
nobody carries together it returned a probability of 1.0 of compound heterozygosity (R16). The 1.0
reflects missing data in the service, so we set its output aside.

### 6.3 Other limits, and what would change the answer

| Limit | Effect | What would address it |
|---|---|---|
| The search window is 170 kb; the gene's regulatory reach is wider (4.3) | A distant regulatory variant remains possible | Check the child's genotypes in the wider regions, or measure *BUB1B* expression |
| Short reads detect substitutions, small indels and exon-scale depth changes (4.3) | A balanced rearrangement, deletion wholly inside an intron, repeat expansion or imprinting defect that explains the condition would displace our answer | Long reads, optical mapping, methylation testing |
| The aneuploidy screen detects gains in about 13% of cells or more | The diagnosis rests on the genotype; the aneuploidy itself is unobserved | Karyotype or single-cell sequencing of cultured cells |
| Mitochondrial calls reach down to about 21% of copies | A low-level mitochondrial cause remains possible | A pipeline that corrects for nuclear copies |
| The data is one sample and one tissue, without RNA, parents or a functional test | Allele 2 and phase stay open. Parental testing showing both variants on one copy with the other copy intact would remove the explanation | Parental testing; BubR1 protein in the child's cultured cells |
| ClinVar and gnomAD are snapshots (ClinVar 2026-09-13) | Classifications can change | Rerun `run_all.sh` before using the result |
| A test showing that p.Asn1002Lys makes normal levels of normal BubR1 | Would leave allele 1 as a carrier state and reopen the search for the second hit | A protein assay |

---

## 7. Secondary finding

The organizers invite secondary findings. Our genome-wide ClinVar sweep (section 4.4) found
six pathogenic or likely pathogenic variants with good read support: allele 1 in *BUB1B* (the
primary finding), the *LZTR1* variant reported below, and four single heterozygous variants in
genes unrelated to the condition.

**`chr22:20996720 C>G`, *LZTR1* p.Tyr748Ter** (transcript `ENST00000646124`). ClinVar 1409252,
Pathogenic/Likely pathogenic, several submitters, no conflict. Heterozygous, with 24 of 48
reads carrying it. It is a stop-gain 162 bases before the last exon-exon junction, so we expect
the cell to destroy the message. We re-derived the change by hand (`consequence-verification.txt`).

* **Possible meaning.** Heterozygous loss of *LZTR1* function predisposes to schwannomatosis,
  a dominant condition of nerve-sheath tumours. Not every carrier develops it: some carrier
  parents in the first report had no symptoms (Piotrowski 2014).
* **Relation to the child's condition.** ClinVar also links this variant to Noonan syndrome,
  which, like this child's condition, comes with a tumour predisposition. Noonan syndrome 10
  has been described with dominant missense variants in *LZTR1* (Yamamoto 2015). Noonan
  syndrome 2 needs damaging variants in both copies (Johnston 2018), and the child has one
  damaging *LZTR1* allele. We attribute the child's condition to the *BUB1B* genotype; whether
  this variant modifies it is unexamined.
* **Next step.** This is a research observation. Confirmation in an accredited laboratory and
  discussion with the treating team and a genetic counsellor come before any action.

We report the variants that bear on the child's presentation, and ClinVar links this one to a
condition that overlaps with it. Sharing health information about the four unrelated variants is
the family's decision.

---

## 8. Additional findings

These findings are extra to the answer. Each may help with the next case or with Track 2.

### 8.1 AlphaMissense and the missense variant

AlphaMissense (Cheng 2023) has scored every possible amino-acid change in every human
protein, and the tables are public. For p.Asn1002Lys the score is **0.9229**, labelled
"likely pathogenic" by its developers. That is the 96th percentile of the 6,940 changes in
*BUB1B* reachable by one nucleotide change (`alphamissense-bub1b.txt`). The score supports a
damaging effect, and three results qualify it.

* The score comes mostly from the position. Eighteen of the 19 possible changes at residue
  1002 score above the developers' cutoff of 0.564, and lysine is fifth. So most of the 0.92
  is the residue.
* A high score says little about how much protein is left. Of five *BUB1B* missense alleles
  whose protein is low in the laboratory (Suijkerbuijk 2010), AlphaMissense scores four below its
  cutoff (p.Arg727Cys 0.24, p.Arg814His 0.22, p.Leu844Phe 0.29, p.Ile909Thr 0.45) and one
  above (p.Leu1012Pro 0.87). A low score would have left our variant possible.
* Under ClinGen's calibrated bands (Bergquist 2025) 0.9229 is "moderate", worth +2 points.
  Allele 2 would have 3 points, still a VUS (section 5.3).

We treat the score as supporting context and claim no ACMG criterion from it. In our Track 2
work (R20), twelve protein language models put the change between the 14th and 79th percentile of
damage for this gene and rated lysine among the mildest changes at residue 1002, and four
structure predictors returned the unchanged fold. None separated known disease alleles from
neutral ones, so we use none as evidence.

### 8.2 Allele-balance screens have little power for aneuploidy in this disease

Allele-balance and depth screens look for a chromosome that is gained or lost in a large share
of cells. In MVA different cells carry different errors. If one cell in three carried an error,
spread over 23 chromosomes, each chromosome would be wrong in about 1.4% of cells, far below the
13% this test can see. A negative result is therefore uninformative about the disease, and the
same arithmetic applies to a variance statistic. Allele-balance spread was 1.27 times the
binomial expectation, tracked chromosome GC content, and had one outlier confined to part of
chromosome 20, which points to a technical origin (R19). Neither the per-chromosome shift nor
the variance statistic can detect a mosaic at this level.

### 8.3 The mitochondrial genome and the 2% cutoff

The supplied call set lacks the mitochondrial genome, which the raw reads hold. We streamed 273
million reads, kept the 0.07% that match the 16,569-base mitochondrial genome, and covered every
position at a median of 1,609 times. A fixed 2% cutoff for low-level variants listed 1,034
"variants", two of which matched pathogenic records (R21). They were sequencing error and
copies of mitochondrial sequence in the nucleus. At this coverage a fixed cutoff is arbitrary, so
we set it from the measured noise: the 99.9th percentile of the background, 21.2%, where a 2%
cutoff would admit 555 positions (`mitochondrial-screen.txt`). The variants that remain are the
common ones that define a mitochondrial lineage, and six universal positions served as a positive
control. The method calls variants down to about 21%.

### 8.4 ClinVar classifies few of this gene's missense variants

Of the 1,490 *BUB1B* missense records in ClinVar, 1,426 (96%) are of uncertain significance,
35 are benign or likely benign, 28 conflicting and one likely pathogenic (R12). A well-studied
MVA1 missense variant is itself a VUS. Partly working alleles need laboratory data to be
classed, so they accumulate as VUS. To estimate how reliable predictors are for these genes we
pooled 24 spindle-checkpoint genes: 14 pathogenic or likely pathogenic missense variants against
300 benign or likely benign ones, all scored with AlphaMissense. The estimate is too wide to use,
with a bootstrap interval on the likelihood ratio from 0.0 to 35.7 (`family-calibration.txt`).
The limit is the number of labelled examples. We built this calibration after seeing the
variant's score. ClinGen requires the tool and threshold to be chosen in advance, so this
calibration serves as context only.

### 8.5 The genotype and Track 2

* Allele 1 gives no protein. A treatment that works by making ribosomes read through a stop
  signal needs a message to act on, and the stop is 747 bases before the last junction, so
  little message should be left.
* Allele 2 makes protein, probably in reduced amount, and nobody has measured how much.
  Raising the level or the stability of that protein is the direction the genotype points to.
* How much protein is enough is unknown. The titration in section 5.2 comes from one cultured
  cell line, and the BubR1 level at which a person develops disease is unmeasured.

These points follow from the genotype and serve as starting hypotheses for Track 2, whose
submission covers candidate medicines and carries its own disclaimer.

---

## 9. Reproducibility, scalability, privacy and use of AI

**Code.** The repository holds the 24 scripts, 5 package modules and 7 test files that produced
the results. `MANIFEST.md` lists each file and what it does. The analysis
scripts use the Python standard library plus `bcftools`, `minimap2`, `samtools` and `bbmap`.
`./run_all.sh` regenerates the results on a machine with the data. `make check` runs the
offline tests.

**Hardware and time.** The main path from the VCF to the predictions file ran on 4 vCPUs, 16
GB of memory and 30 GB of disk, with no GPU, in about 55 minutes on an earlier version of the pipeline (R17),
most of it downloading public references. The raw-read checks need the 85 GB of FASTQ files.

**Results in the repository.** The files in `results/` cover most numbers in this report. R1 to
R21 are listed in `results/RECORDED-CHECKS.md` (section 3).

**Scalability.** The ranking takes no gene, disease or phenotype as input (it uses ClinVar
labels). The annotation step (`bcftools csq`) holds one chromosome's reference at a time, so
memory and disk stay constant as the number of genomes grows. The same scripts run on any
single-sample GRCh38 VCF.

**Privacy.** Variant-effect annotation, ClinVar, gnomAD and the predictors ran against bulk
public files that we downloaded in full and queried on the local machine. The gnomAD browser and
Ensembl received the identifiers of the two reported variants and the secondary finding; the
structure and language-model services named under "Use of generative AI" received public
identifiers and a small number of candidate sites. Results that describe the child pass through `scripts/public_view.py`, which
withholds other genes' genotypes, rsIDs, feature names and mitochondrial positions. The repository holds
code and aggregate results. All copies of the data will be deleted within 30 days of the end of
the competition.

**Use of generative AI.** We used Anthropic Claude (Claude Code) throughout, under commercial
terms with no training on customer content, for analysis design, code, interpretation and
drafting.

---

## 10. Tools, data and references

### Tools and data

Links to every tool, dataset and paper are in `docs/references.md` in the public repository.

| Used for | Name and version | Reference |
|---|---|---|
| Effect of each variant on the protein | `bcftools csq` 1.19 to 1.24 with Ensembl release 116 | Danecek & McCarthy 2017, *Bioinformatics* 33:2037 |
| Re-alignment of raw reads | `minimap2`, `bbmap` 40.02, `mappy` 2.31 | Li 2018, *Bioinformatics* 34:3094 |
| Splice prediction | SpliceAI 1.3.1, run locally | Jaganathan 2019, *Cell* 176:535 |
| Variant frequencies | gnomAD v4.1 | Chen 2024, *Nature* 625:92 |
| Variant classifications | ClinVar, snapshot 2026-09-13 | ncbi.nlm.nih.gov/clinvar |
| Gene–disease rating and dosage curation | ClinGen, *BUB1B* | search.clinicalgenome.org/kb/genes/HGNC:1149 and /kb/gene-dosage/HGNC:1149 |
| Phenotype annotations | HPO release 2026-09-02 | Köhler 2021, *Nucleic Acids Res* 49:D1207 |
| Missense score: AlphaMissense | AlphaMissense | Cheng 2023, *Science* 381:eadg7492 |
| Missense score: REVEL | REVEL v1.3 | Ioannidis 2016, *Am J Hum Genet* 99:877 |
| Missense score: VEST4 | VEST4, through the UCSC dbNSFP track | Carter 2013, *BMC Genomics* 14 (Suppl 3):S3 |
| Missense score: PolyPhen-2 | PolyPhen-2, through the UCSC dbNSFP track | Adzhubei 2010, *Nat Methods* 7:248 |
| Missense score: SIFT | SIFT, through the UCSC dbNSFP track | Ng & Henikoff 2003, *Nucleic Acids Res* 31:3812 |
| Protein structure | AlphaFold model AF-O60566-F1 v6; UniProt O60566; Biopython | Jumper 2021, *Nature* 596:583 |

Also used: Python 3.11, `samtools`, GTEx (for the reach of the gene's regulatory signals).

### Papers our reasoning rests on

Ten of the 13 papers are free in PubMed Central (2 to 10 and 12); three are paywalled (1, 11 and
13). `docs/references.md` in the public repository gives the access status of each.

1. Hanks S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic
   mutations in *BUB1B*. *Nat Genet* 2004;36:1159. PMID 15475955.
2. Suijkerbuijk SJE, et al. Molecular causes for BUBR1 dysfunction in the human cancer
   predisposition syndrome mosaic variegated aneuploidy. *Cancer Res* 2010;70:4891.
   PMID 20516114.
3. Ochiai H, et al. TALEN-mediated single-base-pair editing identification of an intergenic
   mutation upstream of *BUB1B* as causative of PCS (MVA) syndrome. *Proc Natl Acad Sci USA*
   2014;111:1461. PMID 24344301. We cite it as the reason a regulatory second hit is a real
   possibility.
4. Richards S, et al. Standards and guidelines for the interpretation of sequence variants.
   *Genet Med* 2015;17:405. PMID 25741868.
5. Tavtigian SV, et al. Fitting a naturally scaled point system to the ACMG/AMP variant
   classification guidelines. *Hum Mutat* 2020;41:1734. PMID 32720330.
6. Abou Tayoun AN, et al. Recommendations for interpreting the loss of function PVS1
   ACMG/AMP variant criterion. *Hum Mutat* 2018;39:1517. PMID 30192042.
7. Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity
   classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet*
   2022;109:2163. PMID 36413997.
8. Bergquist T, et al. Calibration of additional computational tools expands ClinGen
   recommendation options for variant classification with PP3/BP4 criteria. *Genet Med*
   2025;27:101402. PMID 40084623.
9. Brnich SE, et al. Recommendations for application of the functional evidence PS3/BS3
   criterion. *Genome Med* 2019;12:3. PMID 31892348.
10. Piotrowski A, et al. Germline loss-of-function mutations in *LZTR1* predispose to an
    inherited disorder of multiple schwannomas. *Nat Genet* 2014;46:182. PMID 24362817.
11. Yamamoto GL, et al. Rare variants in *SOS2* and *LZTR1* are associated with Noonan
    syndrome. *J Med Genet* 2015;52:413. PMID 25795793.
12. Johnston JJ, et al. Autosomal recessive Noonan syndrome associated with biallelic *LZTR1*
    variants. *Genet Med* 2018;20:1175. PMID 29469822.
13. Bohers E, et al. Gradual reduction of BUBR1 protein levels results in premature
    sister-chromatid separation then in aneuploidy. *Hum Genet* 2008;124:473. PMID 18932004.
