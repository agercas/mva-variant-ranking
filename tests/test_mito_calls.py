"""The sentence after the seam artefacts must follow from the calls, not be written in advance.

Calls are ``(position, ref, alt, fraction, depth)`` tuples. The ones below are synthetic.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_script():
    spec = importlib.util.spec_from_file_location(
        "mito_call_variants", ROOT / "scripts" / "mito_call_variants.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_sub_homoplasmic_calls_are_counted() -> None:
    rest = [(100 + i, "A", "G", 0.99, 500) for i in range(10)] + [
        (300, "C", "T", 0.40, 500),
        (301, "C", "T", 0.60, 500),
    ]
    assert load_script().remaining_calls_sentence(rest) == (
        "12 calls remain: 10 homoplasmic and 2 sub-homoplasmic (below 95%)."
    )


def test_the_homoplasmic_claim_is_made_only_when_it_is_true() -> None:
    rest = [(100 + i, "A", "G", 0.97, 500) for i in range(3)]
    assert load_script().remaining_calls_sentence(rest) == "every remaining call is homoplasmic: 3 variants."
    assert "every remaining call" not in load_script().remaining_calls_sentence(
        [*rest, (400, "T", "C", 0.949, 500)]
    )
