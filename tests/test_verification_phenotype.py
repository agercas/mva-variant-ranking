"""The reverse phenotype check must run from a terms file and must not carry the supplied terms.

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
name: Term A

[Term]
id: HP:0000101
name: Term B
is_a: HP:0000100 ! Term A

[Term]
id: HP:0000200
name: Term C
is_a: HP:0000001 ! All

[Term]
id: HP:0000201
name: Term D
is_a: HP:0000200 ! Term C
"""

HPOA = "\n".join(
    [
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect",
        "OMIM:257300\tMosaic variegated aneuploidy syndrome 1\t\tHP:0000201\tX\tTAS\t\t\t\t\tP",
        "OMIM:257300\tMosaic variegated aneuploidy syndrome 1\t\tHP:0000101\tX\tTAS\t\t\t\t\tP",
        "OMIM:900001\tSynthetic disease one\t\tHP:0000201\tX\tTAS\t\t\t\t\tP",
        "OMIM:900002\tSynthetic disease two\t\tHP:0000100\tX\tTAS\t\t\t\t\tP",
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
    assert "Term D" in out.stdout  # labels are read from the ontology, not from the source


def test_the_public_view_accepts_the_scripts_own_output(tmp_path: Path) -> None:
    """The view anchors on the table header the script prints; a renamed header must not slip past."""
    raw = run(tmp_path, ["HP:0000201", "HP:0000101"])
    assert raw.returncode == 0, raw.stderr
    view = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "public_view.py"), "--kind", "phenotype"],
        input=raw.stdout,
        capture_output=True,
        text=True,
    )
    assert view.returncode == 0, view.stderr
    assert "Term D" not in view.stdout and "Term B" not in view.stdout
    assert "term names are not shown; counts only" in view.stdout
    assert "Link grades" in view.stdout


def test_a_missing_second_terms_file_is_said_in_the_output(tmp_path: Path) -> None:
    """run_all.sh discards standard error, so the note must be on standard output."""
    out = run(tmp_path, ["HP:0000201"])
    assert out.returncode == 0, out.stderr
    assert "NOTE: no second terms file was found; running on the 1 terms in terms.json only." in out.stdout
    assert str(tmp_path) not in out.stdout  # no machine path in a file that may be published


def load_script():
    import importlib.util

    spec = importlib.util.spec_from_file_location("phenotype_reverse_check", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_the_shipped_watch_list_names_only_the_mva_related_diseases() -> None:
    script = load_script()
    assert all("MVA" in label or "Atelis" in label for label in script.WATCH.values())
    assert len(script.WATCH) == 6


def test_further_diseases_come_from_the_private_file_after_the_shipped_ones(tmp_path: Path) -> None:
    script = load_script()
    assert script.load_watch(tmp_path / "missing.json") == script.WATCH
    extra = tmp_path / "extra.json"
    extra.write_text('{"OMIM:900001": "Synthetic disease one", "OMIM:257300": "an older label"}')
    watch = script.load_watch(extra)
    assert list(watch)[: len(script.WATCH)] == list(script.WATCH)
    assert watch["OMIM:900001"] == "Synthetic disease one"
    assert watch["OMIM:257300"] == script.WATCH["OMIM:257300"]  # the shipped label wins
    extra.write_text("[]")
    try:
        script.load_watch(extra)
    except ValueError as err:
        assert "JSON object" in str(err)
    else:
        raise AssertionError("a file that is not an object must be refused")


def test_stops_clearly_without_a_terms_file(tmp_path: Path) -> None:
    out = run(tmp_path, None)
    assert out.returncode == 1
    assert "no HPO terms found" in out.stdout


def test_source_does_not_hard_code_the_supplied_terms() -> None:
    text = SCRIPT.read_text()
    assert set(re.findall(r"HP:\d{7}", text)) <= {"HP:0000001", "HP:0000002"}  # placeholders only
