"""Contig conventions are the highest-risk detail in Track 1 — test them hard."""

from __future__ import annotations

import pytest

from mva.genomics.coordinates import (
    CoordinateError,
    Variant,
    normalise_contig,
    to_submission_contig,
    to_vcf_contig,
)


@pytest.mark.parametrize(
    ("given", "expected"),
    [("chr15", "15"), ("15", "15"), ("CHR15", "15"), ("chrX", "X"), ("x", "X"), ("chrMT", "MT")],
)
def test_normalise_contig_accepts_every_convention(given: str, expected: str) -> None:
    assert normalise_contig(given) == expected


@pytest.mark.parametrize("bad", ["chr23", "99", "scaffold_1", "", "chrZ"])
def test_normalise_contig_rejects_non_canonical(bad: str) -> None:
    with pytest.raises(CoordinateError):
        normalise_contig(bad)


def test_submission_form_always_has_chr_prefix() -> None:
    # The whole point: the VCF is bare, the submission template is prefixed.
    assert to_vcf_contig("chr15") == "15"
    assert to_submission_contig("15") == "chr15"
    assert to_submission_contig("chr15") == "chr15"


def test_variant_equality_is_convention_independent() -> None:
    assert Variant("chr15", 40163000, "C", "T") == Variant("15", 40163000, "C", "T")
    assert len({Variant("chr15", 1, "A", "G"), Variant("15", 1, "a", "g")}) == 1


@pytest.mark.parametrize("text", ["chr15-40163000-C-T", "15:40163000:C:T", "15-40163000-c-t"])
def test_variant_parse_accepts_both_separators_and_case(text: str) -> None:
    assert Variant.parse(text) == Variant("15", 40163000, "C", "T")


def test_variant_rejects_symbolic_alleles() -> None:
    # The Track 1 schema has no representation for these, so fail loudly.
    with pytest.raises(CoordinateError):
        Variant("15", 100, "A", "<DEL>")


def test_variant_rejects_ref_equal_alt() -> None:
    with pytest.raises(CoordinateError):
        Variant("15", 100, "A", "A")


def test_variant_rejects_non_positive_position() -> None:
    with pytest.raises(CoordinateError):
        Variant("15", 0, "A", "G")


def test_external_id_formats(variant_a: Variant) -> None:
    assert variant_a.to_gnomad_id() == "15-40163000-C-T"  # gnomAD: no chr
    assert variant_a.to_myvariant_id() == "chr15:g.40163000C>T"  # MyVariant: chr
    assert variant_a.submission_contig == "chr15"


def test_indel_myvariant_id_uses_delins() -> None:
    indel = Variant("15", 100, "AT", "G")
    assert indel.to_myvariant_id() == "chr15:g.100_101delinsG"
    assert not indel.is_snv
