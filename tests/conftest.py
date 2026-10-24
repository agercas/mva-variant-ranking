"""Shared pytest configuration.

🔴 No fixture in this repository may contain real proband data. Everything here
is synthetic, and the coordinates are deliberately arbitrary.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from mva.genomics.coordinates import Variant  # noqa: E402
from mva.genomics.submission import PROBAND_ID, Prediction  # noqa: E402


@pytest.fixture
def variant_a() -> Variant:
    """A synthetic variant. Not from any real dataset."""
    return Variant.parse("15-40163000-C-T")


@pytest.fixture
def variant_b() -> Variant:
    """A second synthetic variant, for compound-het tests."""
    return Variant.parse("15-40190000-G-A")


@pytest.fixture
def variant_c() -> Variant:
    """A third synthetic variant, for tests that need a wrong pairing partner."""
    return Variant.parse("15-40200000-T-A")


@pytest.fixture
def pair_prediction(variant_a: Variant, variant_b: Variant) -> Prediction:
    """A compound-heterozygous hypothesis expressed as a single row."""
    return Prediction(
        proband_id=PROBAND_ID,
        variant_1=variant_a,
        variant_2=variant_b,
        epcr=0.95,
        finding_type="primary",
        notes="synthetic test fixture",
    )
