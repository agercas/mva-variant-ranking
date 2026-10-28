"""The competing-diagnosis list must not grow inside the shipped source.

Further genes are read from a git-ignored file at run time. These tests use synthetic gene
symbols only.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def scan():
    return load_script("differential_and_secondary_scan")


def test_absent_extra_file_adds_nothing(scan, tmp_path: Path) -> None:
    assert scan.load_differential(tmp_path / "missing.json") == scan.DIFFERENTIAL


def test_extra_genes_come_first_so_the_output_order_is_unchanged(scan, tmp_path: Path) -> None:
    extra = tmp_path / "extra.json"
    extra.write_text(json.dumps({"GENEA": "Some synthetic syndrome", "GENEB": "Another"}))
    genes = scan.load_differential(extra)
    assert list(genes)[:2] == ["GENEA", "GENEB"]
    assert list(genes)[2:] == list(scan.DIFFERENTIAL)
    assert "GENEA" not in scan.DIFFERENTIAL


def test_the_shipped_description_wins_over_the_private_one(scan, tmp_path: Path) -> None:
    gene = next(iter(scan.DIFFERENTIAL))
    extra = tmp_path / "extra.json"
    extra.write_text(json.dumps({gene: "an older description"}))
    assert scan.load_differential(extra)[gene] == scan.DIFFERENTIAL[gene]


def test_only_the_mva_related_and_reported_genes_are_named(scan) -> None:
    """The competing syndromes other than the MVA-related genes are counted, not named.

    Adding a gene here, or to the public view's whitelist, is a decision about what the
    release says about the child; this test makes it a deliberate one.
    """
    view = load_script("public_view")
    assert set(scan.DIFFERENTIAL) == {"LZTR1", "CEP57", "TRIP13", "BUB1", "BUB1B"}
    assert {"CEP57", "TRIP13", "BUB1"} == view.NAMED_DIFFERENTIALS
    assert {"BUB1B", "LZTR1"} == view.REPORTED_GENES
    assert view.NAMED_DISEASES == ("Mosaic variegated", "MVA", "Atelis")


def test_withheld_row_markers_do_not_pile_up() -> None:
    """Hiding rows of a table that already has a marker must leave one marker, and a second pass changes nothing."""
    view = load_script("public_view")
    text = (
        "   rank    score  disease\n"
        "      1     56.9  Synthetic Rival Syndrome\n"
        "  (other diseases withheld)\n"
        "      3     54.8  Mosaic variegated aneuploidy syndrome 1\n"
    )
    once = view.keep_named_disease_rows(text)
    assert once.count("(other diseases withheld)") == 1
    assert "Synthetic Rival" not in once and "Mosaic variegated" in once
    assert view.keep_named_disease_rows(once) == once


def test_a_file_that_is_not_an_object_is_refused(scan, tmp_path: Path) -> None:
    extra = tmp_path / "extra.json"
    extra.write_text('["GENEA"]')
    with pytest.raises(ValueError, match="JSON object"):
        scan.load_differential(extra)


def test_the_output_says_when_no_extra_gene_was_loaded(scan, tmp_path: Path) -> None:
    missing = scan.load_differential(tmp_path / "missing.json")
    assert "no extra differential genes were loaded" in scan.no_extra_genes_note(missing)
    empty = tmp_path / "empty.json"
    empty.write_text("{}")
    assert scan.no_extra_genes_note(scan.load_differential(empty)) is not None
    extra = tmp_path / "extra.json"
    extra.write_text('{"GENEA": "x"}')
    assert scan.no_extra_genes_note(scan.load_differential(extra)) is None


def test_the_shipped_list_names_only_what_the_public_view_names(scan) -> None:
    """Anything in the script but not in public_view's named set would defeat its redaction."""
    view = load_script("public_view")
    assert set(scan.DIFFERENTIAL) <= set(view.NAMED_DIFFERENTIALS) | set(view.REPORTED_GENES)
