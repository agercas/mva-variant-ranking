"""The ClinGen-calibrated bands must give the same strength for benign and damaging scores.

Bounds are those of Pejaver et al. 2022, Table 2. Scores are arbitrary points inside each band.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def calibrated():
    spec = importlib.util.spec_from_file_location("insilico_panel", ROOT / "scripts" / "insilico_panel.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules["insilico_panel"] = module  # dataclasses look the module up by name
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop("insilico_panel", None)
    return module.CALIBRATED


@pytest.mark.parametrize(
    ("tool", "score", "expected"),
    [
        ("REVEL", 0.472, "indeterminate"),
        ("REVEL", 0.25, "BP4_supporting"),
        ("REVEL", 0.132, "BP4_moderate"),
        ("REVEL", 0.010, "BP4_strong"),
        ("REVEL", 0.002, "BP4_verystrong"),
        ("REVEL", 0.70, "PP3_supporting"),
        ("REVEL", 0.95, "PP3_strong"),
        ("VEST4", 0.433, "BP4_supporting"),
        ("VEST4", 0.231, "BP4_moderate"),
        ("VEST4", 0.841, "PP3_supporting"),
        ("VEST4", 0.885, "PP3_moderate"),
        ("PolyPhen2", 0.997, "PP3_supporting"),
        ("PolyPhen2", 0.113, "BP4_supporting"),
        ("PolyPhen2", 0.005, "BP4_moderate"),
        ("SIFT", 0.0, "PP3_moderate"),
        ("SIFT", 0.001, "PP3_supporting"),
        ("SIFT", 0.5, "BP4_moderate"),
        ("SIFT", 0.327, "BP4_moderate"),
        ("SIFT", 0.2, "BP4_supporting"),
        ("SIFT", 0.080, "BP4_supporting"),
        ("SIFT", 0.079, "indeterminate"),
        ("SIFT", 0.004, "indeterminate"),
    ],
)
def test_strength_is_read_strongest_first(calibrated, tool: str, score: float, expected: str) -> None:
    assert calibrated[tool].strength(score) == expected


def test_a_missing_score_has_no_strength(calibrated) -> None:
    assert calibrated["REVEL"].strength(None) == "no score"
