"""The reverse phenotype check must run from a terms file and must not carry the child's terms.

The ontology and annotations below are tiny synthetic stand-ins for HPO's ``hp.obo`` and
``phenotype.hpoa``. They contain nothing from the child.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "verification" / "phenotype_reverse_check.py"

OBO = """\
[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000100
name: Neoplasm

[Term]
id: HP:0000101
name: Neoplasm of muscle
is_a: HP:0000100 ! Neoplasm

[Term]
id: HP:0000200
name: Growth abnormality
is_a: HP:0000001 ! All

[Term]
id: HP:0000201
name: Short stature
is_a: HP:0000200 ! Growth abnormality
"""

HPOA = "\n".join(
    [
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect",
        "OMIM:257300\tMosaic variegated aneuploidy syndrome 1\t\tHP:0000201\tX\tTAS\t\t\t\t\tP",
        "OMIM:257300\tMosaic variegated aneuploidy syndrome 1\t\tHP:0000101\tX\tTAS\t\t\t\t\tP",
        "OMIM:218040\tCostello syndrome\t\tHP:0000201\tX\tTAS\t\t\t\t\tP",
        "OMIM:151623\tLi-Fraumeni syndrome\t\tHP:0000100\tX\tTAS\t\t\t\t\tP",
    ]
)


def run(tmp_path: Path, terms: list[str] | None) -> subprocess.CompletedProcess[str]:
    hpo = tmp_path / "hpo"
    hpo.mkdir()
    (hpo / "hp.obo").write_text(OBO)
    (hpo / "phenotype.hpoa").write_text(HPOA)
    args = [
        sys.executable,
        str(SCRIPT),
        "--hpo-dir",
        str(hpo),
        "--withheld-terms",
        str(tmp_path / "none.json"),
    ]
    if terms is not None:
        terms_file = tmp_path / "terms.json"
        terms_file.write_text('["' + '", "'.join(terms) + '"]')
        args += ["--terms", str(terms_file)]
    else:
        args += ["--terms", str(tmp_path / "missing.json")]
    return subprocess.run(args, capture_output=True, text=True)


def test_runs_from_a_terms_file(tmp_path: Path) -> None:
    out = run(tmp_path, ["HP:0000201", "HP:0000101"])
    assert out.returncode == 0, out.stderr
    assert "QUESTION 1" in out.stdout and "QUESTION 2" in out.stdout
    assert "MVA1 (BUB1B)" in out.stdout
    assert "Short stature" in out.stdout  # labels are read from the ontology, not from the source


def test_stops_clearly_without_a_terms_file(tmp_path: Path) -> None:
    out = run(tmp_path, None)
    assert out.returncode == 1
    assert "no HPO terms found" in out.stdout


def test_source_does_not_hard_code_the_childs_terms() -> None:
    text = SCRIPT.read_text()
    assert set(re.findall(r"HP:\d{7}", text)) <= {"HP:0000001", "HP:0000002"}  # placeholders only
