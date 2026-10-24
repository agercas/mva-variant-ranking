"""Offline reimplementation of the Track 1 metrics.

Lets us compare candidate submissions **before** spending one of six uploads.

Unlike the first version of this module, the rules below are **not inferred**.
The organizers publish their scorer in the hackathon Space
(``evaluation.py``, fetched 2026-09-13); this module mirrors its semantics
exactly, including the parts that are surprising. See
the organizers' published rules for the line-by-line comparison.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from mva.genomics.coordinates import Variant
from mva.genomics.submission import PROBAND_ID, Prediction

__all__ = ["PROBAND_ID", "RANK_POINT_TIERS", "Score", "f_max", "rank_points", "score", "tier_points"]

#: Points by rank, as a **step function** — verbatim from the organizers'
#: ``evaluation.py::RANK_POINT_TIERS``. Rank 1 pays 100; ranks 2–3 pay 50;
#: ranks 4–5 pay 25; ranks 6–10 pay 10; rank 11+ pays nothing.
#:
#: 🔴 There is no smooth decay. Rank 2 is worth exactly as much as rank 3, and
#: half of rank 1. Moving a candidate from rank 4 to rank 3 doubles its payout;
#: moving it from rank 2 to rank 3 changes nothing.
RANK_POINT_TIERS: tuple[tuple[int, float], ...] = ((1, 100.0), (3, 50.0), (5, 25.0), (10, 10.0))


def tier_points(rank: int) -> float:
    """Points paid at a 1-based rank.

    Args:
        rank: 1-based position after sorting by ``epcr`` descending.

    Returns:
        The tier payout, or ``0.0`` beyond rank 10.
    """
    for max_rank, points in RANK_POINT_TIERS:
        if rank <= max_rank:
            return points
    return 0.0


@dataclass(frozen=True)
class Score:
    """The outcome of scoring one submission."""

    rank_points: float
    f_max: float
    f_max_threshold: float | None
    full_match_rank: int | None
    partial_match_rank: int | None
    matched: tuple[Variant, ...]
    missed: tuple[Variant, ...]

    @property
    def best_rank(self) -> int | None:
        """The rank that actually paid, full match preferred."""
        return self.full_match_rank if self.full_match_rank is not None else self.partial_match_rank

    @property
    def is_perfect(self) -> bool:
        """True for the saturated leaderboard result: 100.0 / 1.000."""
        return self.rank_points >= 100.0 and self.f_max >= 1.0

    def render(self) -> str:
        """A short human-readable summary."""
        rank = self.best_rank if self.best_rank is not None else "—"
        kind = "full" if self.full_match_rank is not None else "partial"
        threshold = "—" if self.f_max_threshold is None else f"{self.f_max_threshold:.3f}"
        return (
            f"rank points {self.rank_points:6.1f}   F-max {self.f_max:5.3f} "
            f"(at epcr >= {threshold})   best rank {rank} ({kind})   "
            f"matched {len(self.matched)}/{len(self.matched) + len(self.missed)}"
        )


def _ranked(predictions: Sequence[Prediction]) -> list[Prediction]:
    """Sort by epcr descending, ties broken by original file order.

    Mirrors ``sorted(enumerate(rows), key=lambda x: (-x[1][1], x[0]))``.
    """
    return [p for _, p in sorted(enumerate(predictions), key=lambda x: (-float(x[1].epcr), x[0]))]


def _variants_of(prediction: Prediction) -> frozenset[Variant]:
    if prediction.variant_2 is not None:
        return frozenset({prediction.variant_1, prediction.variant_2})
    return frozenset({prediction.variant_1})


def rank_points(
    predictions: Sequence[Prediction], truth: Sequence[Variant]
) -> tuple[float, int | None, int | None]:
    """Compute rank points the way the organizers' scorer does.

    🔴 A **full match requires exact set equality** with the truth set. A row
    holding one correct variant paired with an incorrect second variant is a
    *partial* match, not a full one — and two correct variants split across two
    separate rows are never a full match at all.

    🔴 A full match anywhere in the top 10 **suppresses partial credit
    entirely**, even when the partial would have paid more. A full match at
    rank 5 scores 25; a partial at rank 1 that it outranks would have scored 50.
    This is the organizers' behaviour, not a bug here — do not "fix" it.

    Args:
        predictions: Candidate rows; ranked internally by ``epcr``.
        truth: The causal variant(s). Two entries means compound heterozygous.

    Returns:
        ``(points, full_match_rank, partial_match_rank)``, ranks 1-based or
        ``None``.
    """
    truth_set = frozenset(truth)
    if not truth_set:
        raise ValueError("truth set is empty")

    ordered = _ranked(predictions)
    is_compound_het = len(truth_set) == 2

    full_match_rank: int | None = None
    for rank, prediction in enumerate(ordered, start=1):
        if _variants_of(prediction) == truth_set:
            full_match_rank = rank
            break

    partial_match_rank: int | None = None
    if is_compound_het and full_match_rank is None:
        for rank, prediction in enumerate(ordered, start=1):
            if _variants_of(prediction) & truth_set:
                partial_match_rank = rank
                break

    if full_match_rank is not None:
        return tier_points(full_match_rank), full_match_rank, None
    if partial_match_rank is not None:
        return 0.5 * tier_points(partial_match_rank), None, partial_match_rank
    return 0.0, None, None


def f_max(predictions: Sequence[Prediction], truth: Sequence[Variant]) -> tuple[float, float | None]:
    """Compute F-max and the EPCR threshold that achieves it.

    Swept at the **individual-variant** level: at each distinct threshold the
    call set is the union of variants across every row with ``epcr >=
    threshold``, and F1 is computed against the truth set.

    The practical consequence: rows *below* your best row are free — the sweep
    simply picks a higher threshold that excludes them. Rows *above* it are not:
    they inject false positives into every threshold that admits the correct
    row. So fill all ten slots if you like, but only under your best guess.

    Returns:
        ``(f_max, threshold)``; threshold is ``None`` when nothing ever matched.
    """
    truth_set = frozenset(truth)
    if not truth_set:
        raise ValueError("truth set is empty")

    best_f1, best_threshold = 0.0, None
    for threshold in sorted({float(p.epcr) for p in predictions}, reverse=True):
        called: set[Variant] = set()
        for prediction in predictions:
            if float(prediction.epcr) >= threshold:
                called |= _variants_of(prediction)

        true_positives = len(called & truth_set)
        precision = true_positives / len(called) if called else 0.0
        recall = true_positives / len(truth_set)
        if precision + recall == 0:
            continue
        f1 = 2 * precision * recall / (precision + recall)
        if f1 > best_f1:
            best_f1, best_threshold = f1, threshold

    return best_f1, best_threshold


def score(predictions: Sequence[Prediction], truth: Sequence[Variant]) -> Score:
    """Score a submission against a known truth set.

    🔴 **Secondary findings are scored.** The FAQ says they "won't hurt your
    automated score", but the organizers' code never filters on
    ``finding_type`` — it is carried through as an informational field only.
    A secondary row therefore occupies a rank and contributes variants to the
    F-max sweep like any other. The claim is only true if you keep secondary
    rows *below* every primary one, which is what the official template does
    (0.90 primary, 0.30 secondary) and what ``submission.validate`` enforces.
    """
    truth_set = frozenset(truth)
    points, full_rank, partial_rank = rank_points(predictions, truth)
    best_f1, threshold = f_max(predictions, truth)

    called: set[Variant] = set()
    for prediction in predictions:
        called |= _variants_of(prediction)

    return Score(
        rank_points=points,
        f_max=best_f1,
        f_max_threshold=threshold,
        full_match_rank=full_rank,
        partial_match_rank=partial_rank,
        matched=tuple(sorted(called & truth_set)),
        missed=tuple(sorted(truth_set - called)),
    )
