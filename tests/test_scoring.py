"""The offline scorer must reproduce the organizers' evaluation.py exactly.

Every expectation here was checked against their published implementation by
differential fuzzing (``tests/test_scoring_live.py``), not inferred from the
leaderboard. Where the official behaviour is surprising, the test says so.
"""

from __future__ import annotations

from mva.genomics.coordinates import Variant
from mva.genomics.submission import Prediction
from mva.scoring import f_max, rank_points, score, tier_points


def test_rank_points_are_a_step_function() -> None:
    """Not a smooth decay: 100 / 50 / 50 / 25 / 25 / 10 x5 / 0."""
    assert [tier_points(r) for r in range(1, 12)] == [
        100.0,
        50.0,
        50.0,
        25.0,
        25.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        0.0,
    ]


def test_perfect_paired_call_scores_100(variant_a: Variant, variant_b: Variant) -> None:
    predictions = [Prediction("PROBAND01", variant_a, 0.95, variant_2=variant_b)]
    result = score(predictions, [variant_a, variant_b])
    assert result.rank_points == 100.0
    assert result.f_max == 1.0
    assert result.is_perfect


def test_one_of_a_pair_gets_half_credit(variant_a: Variant, variant_b: Variant) -> None:
    result = score([Prediction("PROBAND01", variant_a, 0.9)], [variant_a, variant_b])
    assert result.rank_points == 50.0
    assert result.partial_match_rank == 1
    assert result.missed == (variant_b,)


def test_a_complete_miss_scores_zero(variant_a: Variant, variant_b: Variant) -> None:
    other = Variant("1", 500, "A", "G")
    result = score([Prediction("PROBAND01", other, 0.9)], [variant_a, variant_b])
    assert result.rank_points == 0.0
    assert result.f_max == 0.0
    assert result.best_rank is None


def test_the_pair_must_share_one_row(variant_a: Variant, variant_b: Variant) -> None:
    """🔴 Split the compound-het pair across two rows and you forfeit half the points.

    A full match requires exact set equality with the truth set, so two correct
    variants on separate rows can only ever register as a partial match — even
    though F-max still reaches 1.0 at the lower threshold.
    """
    split = [
        Prediction("PROBAND01", variant_a, 0.95),
        Prediction("PROBAND01", variant_b, 0.90),
    ]
    paired = [Prediction("PROBAND01", variant_a, 0.95, variant_2=variant_b)]

    split_result = score(split, [variant_a, variant_b])
    assert split_result.rank_points == 50.0
    assert split_result.full_match_rank is None
    assert split_result.f_max == 1.0

    assert score(paired, [variant_a, variant_b]).rank_points == 100.0


def test_a_wrong_partner_downgrades_a_row_to_partial(variant_a: Variant, variant_b: Variant) -> None:
    """One correct variant plus a wrong second is a partial match, not a full one."""
    wrong = Variant("1", 999, "C", "T")
    result = score([Prediction("PROBAND01", variant_a, 0.9, variant_2=wrong)], [variant_a, variant_b])
    assert result.rank_points == 50.0
    assert result.full_match_rank is None


def test_a_full_match_suppresses_a_better_paying_partial(variant_a: Variant, variant_b: Variant) -> None:
    """🔴 Official behaviour, and it is counter-intuitive — do not 'fix' it.

    Once any row matches the pair exactly, partial credit is never consulted,
    even when a partial sitting higher up would have paid more. Here a full
    match at rank 10 pays 10, while the rank-1 partial that would have paid 50
    is ignored. The lesson is not to hedge: put the pair as high as you believe it.
    """
    decoys = [Prediction("PROBAND01", Variant("1", 1000 + i, "A", "G"), 0.9 - 0.05 * i) for i in range(8)]
    predictions = [
        Prediction("PROBAND01", variant_a, 0.99),  # rank 1, would pay 50 as a partial
        *decoys,
        Prediction("PROBAND01", variant_b, 0.10, variant_2=variant_a),  # rank 10, full match
    ]
    result = score(predictions, [variant_a, variant_b])
    assert result.full_match_rank == 10
    assert result.rank_points == 10.0


def test_rows_below_the_answer_are_free(variant_a: Variant, variant_b: Variant) -> None:
    """Padding costs nothing as long as it sits under your best guess.

    F-max sweeps every threshold and keeps the best, so a higher threshold
    simply excludes the padding.
    """
    lean = [Prediction("PROBAND01", variant_a, 0.95, variant_2=variant_b)]
    padded = [
        *lean,
        *[Prediction("PROBAND01", Variant("1", 2000 + i, "A", "G"), 0.5 - 0.01 * i) for i in range(9)],
    ]
    assert score(padded, [variant_a, variant_b]).rank_points == 100.0
    assert score(padded, [variant_a, variant_b]).f_max == score(lean, [variant_a, variant_b]).f_max == 1.0


def test_rows_above_the_answer_are_not_free(variant_a: Variant, variant_b: Variant) -> None:
    """A wrong row ranked above the correct pair injects unavoidable false positives."""
    above = [
        Prediction("PROBAND01", Variant("1", 2000, "A", "G"), 0.99, variant_2=Variant("1", 2001, "C", "T")),
        Prediction("PROBAND01", variant_a, 0.90, variant_2=variant_b),
    ]
    result = score(above, [variant_a, variant_b])
    assert result.rank_points == 50.0  # full match at rank 2
    assert result.f_max < 1.0


def test_ties_at_the_top_depress_f_max(variant_a: Variant, variant_b: Variant) -> None:
    """No threshold can separate rows sharing an epcr, so precision suffers."""
    tied = [
        Prediction("PROBAND01", variant_a, 0.9, variant_2=variant_b),
        Prediction("PROBAND01", Variant("2", 500, "A", "G"), 0.9),
    ]
    separated = [
        Prediction("PROBAND01", variant_a, 0.9, variant_2=variant_b),
        Prediction("PROBAND01", Variant("2", 500, "A", "G"), 0.2),
    ]
    assert f_max(tied, [variant_a, variant_b])[0] < f_max(separated, [variant_a, variant_b])[0]
    assert f_max(separated, [variant_a, variant_b])[0] == 1.0


def test_secondary_findings_are_free_only_below_the_primaries(variant_a: Variant, variant_b: Variant) -> None:
    """🔴 The FAQ says secondary findings cannot hurt you. The code disagrees.

    ``finding_type`` is never read by the organizers' scorer, so a secondary row
    occupies a rank and contributes variants to the F-max sweep like any other.
    It is harmless below the primaries, and costly above them.
    """
    incidental = Variant("13", 32340000, "A", "T")
    base = [Prediction("PROBAND01", variant_a, 0.95, variant_2=variant_b)]

    below = [*base, Prediction("PROBAND01", incidental, 0.05, finding_type="secondary")]
    assert score(below, [variant_a, variant_b]) == score(base, [variant_a, variant_b])

    above = [*base, Prediction("PROBAND01", incidental, 0.99, finding_type="secondary")]
    result = score(above, [variant_a, variant_b])
    assert result.rank_points == 50.0  # the pair is pushed to rank 2
    assert result.f_max < 1.0


def test_rank_points_reports_which_kind_of_match(variant_a: Variant) -> None:
    predictions = [
        Prediction("PROBAND01", Variant("1", 100, "A", "G"), 0.9),
        Prediction("PROBAND01", variant_a, 0.8),
    ]
    points, full_rank, partial_rank = rank_points(predictions, [variant_a])
    assert (points, full_rank, partial_rank) == (50.0, 2, None)


def test_epcr_ties_break_on_file_order(variant_a: Variant, variant_b: Variant) -> None:
    """Equal epcr keeps the submitted order, matching the organizers' stable sort."""
    predictions = [
        Prediction("PROBAND01", Variant("1", 100, "A", "G"), 0.5),
        Prediction("PROBAND01", variant_a, 0.5, variant_2=variant_b),
    ]
    assert score(predictions, [variant_a, variant_b]).full_match_rank == 2


def test_correct_pair_row_below_rank_3_scores_worse_than_omitting_it(
    variant_a: Variant, variant_b: Variant, variant_c: Variant
) -> None:
    """🔴 The trap that makes "fill all ten slots" bad advice.

    Partial credit for recovering one allele of a compound-het pair is a
    *fallback*: ``evaluation.py`` scans for an exact full match, takes the first
    one it finds and stops, and only reaches the 0.5x partial branch when there
    is none. Our top row always contains one true allele whenever the answer key
    shares an allele with it, so it banks 50 points by itself.

    A full match therefore has to *beat* 50 to be worth listing at all, and
    below rank 3 it cannot: rank 4-5 pays 25 and rank 6-10 pays 10. Listing the
    correct pair too low is strictly worse than not listing it.
    """
    truth = [variant_a, variant_c]
    top = Prediction("PROBAND01", variant_a, 0.97, variant_2=variant_b)  # wrong pair, shares variant_a
    correct = Prediction("PROBAND01", variant_a, 0.10, variant_2=variant_c)

    without = score([top, Prediction("PROBAND01", variant_b, 0.5)], truth)
    assert without.rank_points == 50.0, "partial credit on the top row"

    # Same file, plus the CORRECT pair pushed down to rank 4.
    padded = score(
        [
            top,
            Prediction("PROBAND01", variant_b, 0.50),
            Prediction("PROBAND01", variant_c, 0.30),
            correct,
        ],
        truth,
    )
    assert padded.rank_points == 25.0
    assert padded.rank_points < without.rank_points, (
        "adding the correct answer at rank 4 must score worse than omitting it - "
        "if this ever passes, the fallback semantics changed and "
        "the Track 1 rules need revisiting"
    )


def test_correct_pair_row_at_rank_3_is_neutral(
    variant_a: Variant, variant_b: Variant, variant_c: Variant
) -> None:
    """Ranks 2-3 pay the same 50 as the fallback, so a pair row there is safe."""
    truth = [variant_a, variant_c]
    rows = [
        Prediction("PROBAND01", variant_a, 0.97, variant_2=variant_b),
        Prediction("PROBAND01", variant_b, 0.50),
        Prediction("PROBAND01", variant_a, 0.30, variant_2=variant_c),  # correct, rank 3
    ]
    assert score(rows, truth).rank_points == 50.0
