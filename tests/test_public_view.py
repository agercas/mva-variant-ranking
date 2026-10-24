"""The public view of a result file must keep the argument and drop the rest.

Fixtures are small synthetic snippets in the same format as the real outputs. None
of them contains data from the child.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "public_view.py"


def view(kind: str, text: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--kind", kind],
        input=text,
        capture_output=True,
        text=True,
    )


RANKING = """\
 rank  gene           pattern                  alleles  ClinVar P/LP   score
------------------------------------------------------------------------------
    1  BUB1B          het-LoF + het-missense         2             1   115.0
    2  GENEX          hom-LoF                        1             0   100.0
    3  GENEY          hom-LoF                        3             0   100.0

>>> BUB1B: rank 1 of 175  [het-LoF + het-missense, score 115.0]
"""


def test_ranking_names_only_the_reported_gene() -> None:
    out = view("ranking", RANKING)
    assert out.returncode == 0
    assert "BUB1B" in out.stdout
    assert "GENEX" not in out.stdout and "GENEY" not in out.stdout
    assert out.stdout.count("(withheld)") == 2
    assert "hom-LoF" in out.stdout  # the pattern and score stay, so the table is still readable


DIFFERENTIAL = """\
PART 1 - genome-wide ClinVar burden
The proband carries 6 high-quality ClinVar P/LP variants genome-wide.
Heterozygous single P/LP alleles:
  BUB1B      15:40209701 T>G GT=0/1 DP=46  not_provided|Some_disease
  LZTR1      22:20996720 C>G GT=0/1 DP=48  not_provided|Other_disease

PART 2 - differentials
GENEA — Some syndrome - blind rank 1 (kidney, growth)
   no damaging variants at all  →  EXCLUDED at the genotype level

GENEB — Another syndrome
   missense       GT=0/1 rs1111111      ClinVar=Benign   [coords redacted]
   missense       GT=1/1 rs2222222      ClinVar=Benign   [coords redacted]
   all 2 are ClinVar-benign common polymorphisms  →  EXCLUDED

HRAS — Costello syndrome
   no damaging variants at all  →  EXCLUDED at the genotype level

BUB1B — hypothesis
   15:40209701 T>G stop_gained    GT=0/1 DP=46 novel          ClinVar=Pathogenic/Likely_pathogenic
"""


def test_differential_drops_rsids_genotypes_and_unnamed_genes() -> None:
    out = view("differential", DIFFERENTIAL)
    assert out.returncode == 0, out.stderr
    assert "rs1111111" not in out.stdout and "rs2222222" not in out.stdout
    assert "GENEA" not in out.stdout and "GENEB" not in out.stdout
    assert "kidney" not in out.stdout
    assert "2 further named competing diagnoses (genes not named): 1 have no damaging variant" in out.stdout
    assert "1 carry only variants that ClinVar calls benign" in out.stdout
    assert "HRAS — Costello syndrome" in out.stdout
    assert "BUB1B" in out.stdout and "LZTR1" in out.stdout
    assert "DP=46" in out.stdout


MITO = """\
EMPIRICAL NOISE FLOOR — measured, not assumed
  variants above the floor: 3 (2 homoplasmic, 1 sub-homoplasmic)
       263 A>G   98.4%  depth   872
      4769 C>T   21.6%  depth   631  <- circular-origin seam artefact
  5 of 17 calls sit within 25 bp of the circular
  For completeness, 3 P/LP position(s) sit INSIDE the noise
  continuum and are reported as artefacts, not findings:
    15 C>T at 8.0% — below the 21.2% floor
DETECTION LIMIT
"""


def test_mito_view_keeps_the_floor_and_drops_positions() -> None:
    out = view("mito", MITO)
    assert out.returncode == 0, out.stderr
    assert "variants above the floor: 3" in out.stdout
    assert "263 A>G" not in out.stdout and "4769 C>T" not in out.stdout
    assert "15 C>T" not in out.stdout
    assert "positions withheld" in out.stdout


PHENOTYPE = """\
  the child presents with          link        IC  explained via
  Feature one                      DIRECT     1.7  (annotated)
  Feature two                      SPECIFIC   6.2  Some ancestor

  Link grades: 1 direct, 1 specific.
"""


def test_phenotype_view_hides_the_term_names() -> None:
    out = view("phenotype", PHENOTYPE)
    assert out.returncode == 0, out.stderr
    assert "Feature one" not in out.stdout and "Some ancestor" not in out.stdout
    assert "1 direct, 1 specific" in out.stdout


def test_phenotype_view_keeps_only_named_diseases_in_the_ranking() -> None:
    text = (
        PHENOTYPE
        + """
   rank    score  disease
      1     56.9  Costello syndrome
      2     50.0  Some Unnamed Syndrome
      3     54.8  Mosaic variegated aneuploidy syndrome 1   <-- BUB1B
"""
    )
    out = view("phenotype", text)
    assert out.returncode == 0, out.stderr
    assert "Unnamed" not in out.stdout and "Costello" in out.stdout
    assert "(other diseases withheld)" in out.stdout


COMPLETENESS = """\
3. Which inheritance models were actually searched?
  chrY  : 2 damaging variants exist (gene and rsIDs withheld).
          That gene has no OMIM disease association matching this phenotype.
"""


def test_completeness_view_passes_clean_text() -> None:
    out = view("completeness", COMPLETENESS)
    assert out.returncode == 0, out.stderr
    assert "withheld" in out.stdout


def test_completeness_view_refuses_a_stray_rsid() -> None:
    out = view("completeness", COMPLETENESS + "  one is rs123456789\n")
    assert out.returncode == 2 and out.stdout == ""


EXPANDED = """\
ClinVar P/LP records inside the panel spans: 5

BUB1B      chr15    14 variants    4 novel  <-- 1 ClinVar P/LP
GENEZ      chr1      9 variants    0 novel  <-- 1 ClinVar P/LP

=== ClinVar Pathogenic/Likely-pathogenic alleles carried, by gene ===

BUB1B:
   40209701 T>G het DP=46 PASS . Pathogenic

GENEZ:
   123456 A>G HOM DP=30 PASS rs123456789 Pathogenic
"""


def test_expanded_view_shows_alleles_only_for_the_reported_gene() -> None:
    out = view("expanded", EXPANDED)
    assert out.returncode == 0, out.stderr
    assert "40209701 T>G het DP=46" in out.stdout
    assert "rs123456789" not in out.stdout and "123456 A>G" not in out.stdout
    assert "GENEZ: (alleles withheld)" in out.stdout


@pytest.mark.parametrize("kind", ["ranking", "differential", "expanded", "mito", "phenotype", "completeness"])
def test_unexpected_input_is_refused_not_passed_through(kind: str) -> None:
    out = view(kind, "some other text that has none of the anchors\n")
    assert out.returncode == 2
    assert out.stdout == ""
