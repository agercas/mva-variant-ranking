"""Six submissions is not many. Nothing invalid may reach the upload endpoint."""

from __future__ import annotations

from pathlib import Path

import pytest

from mva.genomics.coordinates import Variant
from mva.genomics.submission import (
    MAX_ROWS,
    PROBAND_ID,
    Prediction,
    SubmissionError,
    check_file,
    read,
    suggested_filename,
    validate,
    write,
)


def test_valid_submission_passes(pair_prediction: Prediction) -> None:
    assert validate([pair_prediction]).ok


def test_empty_submission_is_rejected() -> None:
    report = validate([])
    assert not report.ok
    assert "empty" in report.errors[0]


def test_duplicate_epcr_is_an_error(variant_a: Variant, variant_b: Variant) -> None:
    # Ties cannot be separated by any threshold, so they permanently cap F-max.
    predictions = [
        Prediction(PROBAND_ID, variant_a, 0.9),
        Prediction(PROBAND_ID, variant_b, 0.9),
    ]
    report = validate(predictions)
    assert not report.ok
    assert any("duplicate epcr" in error for error in report.errors)


def test_too_many_rows_is_an_error(variant_a: Variant) -> None:
    predictions = [
        Prediction(PROBAND_ID, Variant("1", 1000 + i, "A", "G"), 0.9 - i * 0.05) for i in range(MAX_ROWS + 1)
    ]
    assert any("exceeds" in error for error in validate(predictions).errors)


def test_missing_primary_is_an_error(variant_a: Variant) -> None:
    predictions = [Prediction(PROBAND_ID, variant_a, 0.5, finding_type="secondary")]
    assert any("primary" in error for error in validate(predictions).errors)


def test_multiple_proband_ids_is_an_error(variant_a: Variant, variant_b: Variant) -> None:
    predictions = [Prediction("P1", variant_a, 0.9), Prediction("P2", variant_b, 0.5)]
    assert any("proband_id" in error for error in validate(predictions).errors)


def test_duplicate_hypothesis_is_an_error(variant_a: Variant) -> None:
    predictions = [Prediction(PROBAND_ID, variant_a, 0.9), Prediction(PROBAND_ID, variant_a, 0.5)]
    assert any("duplicate prediction" in error for error in validate(predictions).errors)


def test_secondary_outranking_primary_warns(variant_a: Variant, variant_b: Variant) -> None:
    predictions = [
        Prediction(PROBAND_ID, variant_a, 0.4),
        Prediction(PROBAND_ID, variant_b, 0.8, finding_type="secondary"),
    ]
    report = validate(predictions)
    assert report.ok  # not blocking
    assert any("secondary" in warning for warning in report.warnings)


def test_epcr_out_of_range_is_rejected(variant_a: Variant) -> None:
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(SubmissionError):
            Prediction(PROBAND_ID, variant_a, bad)


def test_pair_with_identical_alleles_is_rejected(variant_a: Variant) -> None:
    with pytest.raises(SubmissionError):
        Prediction(PROBAND_ID, variant_a, 0.9, variant_2=variant_a)


def test_round_trip_writes_chr_prefixed_contigs(tmp_path, pair_prediction: Prediction) -> None:
    path = write([pair_prediction], tmp_path / "sub.csv")
    text = path.read_text()
    assert "chr15" in text
    # Bare contig must not survive into the file; check the data row specifically.
    data_row = text.splitlines()[1]
    assert data_row.split(",")[1] == "chr15"

    restored = read(path)
    assert restored[0].variant_1 == pair_prediction.variant_1
    assert restored[0].variant_2 == pair_prediction.variant_2
    assert check_file(path).ok


def test_check_file_catches_missing_chr_prefix(tmp_path) -> None:
    # Hand-written CSV using the VCF's bare convention — the exact mistake that
    # most plausibly produced the 0.0 scores on the leaderboard.
    path = tmp_path / "bare.csv"
    path.write_text(
        "proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes\n"
        "P,15,40163000,C,T,,,,,0.9,primary,\n"
    )
    report = check_file(path)
    assert not report.ok
    assert any("chr" in error for error in report.errors)


def test_write_refuses_invalid_submission(tmp_path, variant_a: Variant, variant_b: Variant) -> None:
    predictions = [Prediction(PROBAND_ID, variant_a, 0.9), Prediction(PROBAND_ID, variant_b, 0.9)]
    with pytest.raises(SubmissionError):
        write(predictions, tmp_path / "bad.csv")


def test_suggested_filename_matches_organizer_convention() -> None:
    assert suggested_filename("jane-doe", "GATK HaplotypeCaller") == "jane-doe_gatk-haplotypecaller.csv"


def test_a_wrong_proband_id_is_an_error(pair_prediction: Prediction) -> None:
    """The Space accepts only PROBAND01 — the VCF sample name is a different thing.

    Getting this wrong wastes one of six submissions on an 'Unknown proband_id'
    rejection, so it blocks rather than warns.
    """
    from dataclasses import replace

    report = validate([replace(pair_prediction, proband_id="SAMPLE_X")])
    assert not report.ok
    assert any("PROBAND01" in error for error in report.errors)


# ---------------------------------------------------------------------------
# 🔴 Pair-vs-single confidence separation
#
# The check that catches "this pair may be scored as two single guesses" used an
# ABSOLUTE epsilon (|a - b| < 0.1). That is scale-dependent and backwards: once
# epcr values are well calibrated the hedges are small (0.06, 0.04, 0.02)
# and every neighbour falls inside 0.1, so the rule fired on a CORRECT ladder —
# while the same ladder inflated to 0.25/0.13 passed silently. These tests pin
# the scale-free behaviour.
# ---------------------------------------------------------------------------


def _pair_and_single(variant_a: Variant, variant_b: Variant, pair_epcr: float, single_epcr: float):
    return [
        Prediction(PROBAND_ID, variant_a, pair_epcr, "primary", variant_b),
        Prediction(PROBAND_ID, variant_a, single_epcr, "primary"),
    ]


def _separation_warnings(report) -> list[str]:
    return [w for w in report.warnings if "compound-het pair" in w]


def test_warns_when_pair_and_single_are_indistinguishable(variant_a: Variant, variant_b: Variant) -> None:
    """The real failure: a pair and one of its alleles at the same confidence."""
    report = validate(_pair_and_single(variant_a, variant_b, 0.50, 0.49))
    assert _separation_warnings(report), "a 2% separation must still warn"


def test_no_warning_for_small_hedges(variant_a: Variant, variant_b: Variant) -> None:
    """0.06 vs 0.04 is a 1.5x ratio — distinguishable, and the real ladder."""
    report = validate(_pair_and_single(variant_a, variant_b, 0.04, 0.06))
    assert not _separation_warnings(report)
    assert report.ok


def test_separation_check_is_scale_free(variant_a: Variant, variant_b: Variant) -> None:
    """The same ratio must behave the same way at any magnitude."""
    small = _separation_warnings(validate(_pair_and_single(variant_a, variant_b, 0.04, 0.06)))
    large = _separation_warnings(validate(_pair_and_single(variant_a, variant_b, 0.40, 0.60)))
    assert (not small) and (not large)

    tight_small = _separation_warnings(validate(_pair_and_single(variant_a, variant_b, 0.040, 0.044)))
    tight_large = _separation_warnings(validate(_pair_and_single(variant_a, variant_b, 0.400, 0.440)))
    assert tight_small and tight_large


def test_the_shipped_ladder_is_clean(tmp_path: Path) -> None:
    """The CSV we upload must emit no warnings, and must still be what the builder produces.

    What this guards:

    - It builds into a temporary directory, so running the offline suite never
      rewrites the submitted file. ``--out`` is required.
    - ``submission/bub1b-compound-het.csv``, the file that goes into the submission
      form, must equal the builder's output in every column, notes included.

    Comparing ``to_row()`` rather than bytes keeps this insensitive to the CRLF
    that ``csv`` writes and ``.gitattributes`` normalises away on commit.
    """
    import subprocess
    import sys

    rebuilt = tmp_path / "rebuilt.csv"
    subprocess.run(
        [sys.executable, "scripts/build_track1_submission.py", "--out", str(rebuilt)],
        capture_output=True,
        check=True,
    )

    rebuilt_rows = read(rebuilt)
    report = validate(rebuilt_rows)
    assert report.ok, report.errors
    assert not report.warnings, report.warnings

    expected = [prediction.to_row() for prediction in rebuilt_rows]
    root = Path(__file__).resolve().parent.parent

    path = root / "submission/bub1b-compound-het.csv"
    assert path.exists(), "no uploaded predictions CSV found"
    assert [prediction.to_row() for prediction in read(path)] == expected, (
        f"{path} has drifted from scripts/build_track1_submission.py: regenerate it"
    )
