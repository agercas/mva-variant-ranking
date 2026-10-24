"""Track 1 submission file construction and validation.

Six submissions is not many, and two leaderboard entries scored 0.0 — almost
certainly a coordinate-convention mismatch. Nothing is uploaded that has not
passed :func:`validate`.

The rules encoded here come from the organizers' Track 1 rules and the EPCR
policy in the organizers' published rules.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mva.genomics.coordinates import CoordinateError, Variant

#: Column order required by the organizers' template.
COLUMNS: tuple[str, ...] = (
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
)

#: How close, *relative to the larger value*, a pair row and a standalone row of
#: one of its alleles may sit before we warn that a threshold sweep cannot tell
#: them apart. Relative rather than absolute on purpose: an absolute epsilon is
#: scale-dependent, so it fires on well-calibrated small hedges and stays
#: quiet on inflated ones — rewarding overconfidence.
PAIR_SINGLE_RELATIVE_EPSILON: float = 0.25

MAX_ROWS = 10

#: The only proband identifier the scorer accepts. Verified against the Space's
#: ``tabs/submit_track1.py``: "Only PROBAND01 is accepted in this challenge."
#: The VCF sample name (sample name) is a different thing entirely.
PROBAND_ID = "PROBAND01"
FINDING_TYPES = frozenset({"primary", "secondary"})


class SubmissionError(ValueError):
    """Raised when a submission would be invalid."""


@dataclass
class Prediction:
    """One row of the Track 1 predictions file.

    A compound-heterozygous hypothesis is **one row carrying both alleles**, not
    two rows — see the organizers' published rules.
    """

    proband_id: str
    variant_1: Variant
    epcr: float
    finding_type: str = "primary"
    variant_2: Variant | None = None
    notes: str = ""

    def __post_init__(self) -> None:
        if not 0.0 < float(self.epcr) <= 1.0:
            raise SubmissionError(f"epcr must be in (0, 1], got {self.epcr}")
        if self.finding_type not in FINDING_TYPES:
            raise SubmissionError(
                f"finding_type must be one of {sorted(FINDING_TYPES)}, got {self.finding_type!r}"
            )
        if self.variant_2 is not None and self.variant_2 == self.variant_1:
            raise SubmissionError(
                "variant_2 duplicates variant_1; a compound-het pair needs two distinct variants"
            )

    @property
    def is_pair(self) -> bool:
        """True when this row expresses a compound-heterozygous hypothesis."""
        return self.variant_2 is not None

    def to_row(self) -> dict[str, Any]:
        """Render to the template's column layout, with ``chr``-prefixed contigs."""
        row: dict[str, Any] = {
            "proband_id": self.proband_id,
            "chrom_1": self.variant_1.submission_contig,
            "pos_1": self.variant_1.position,
            "ref_1": self.variant_1.ref,
            "alt_1": self.variant_1.alt,
            "chrom_2": "",
            "pos_2": "",
            "ref_2": "",
            "alt_2": "",
            "epcr": f"{float(self.epcr):.6g}",
            "finding_type": self.finding_type,
            "notes": self.notes,
        }
        if self.variant_2 is not None:
            row.update(
                chrom_2=self.variant_2.submission_contig,
                pos_2=self.variant_2.position,
                ref_2=self.variant_2.ref,
                alt_2=self.variant_2.alt,
            )
        return row


@dataclass
class ValidationReport:
    """Outcome of validating a set of predictions."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """True when there are no blocking errors."""
        return not self.errors

    def render(self) -> str:
        """A human-readable report."""
        lines: list[str] = []
        for error in self.errors:
            lines.append(f"  ERROR   {error}")
        for warning in self.warnings:
            lines.append(f"  WARNING {warning}")
        if not lines:
            lines.append("  OK — submission passes every check.")
        return "\n".join(lines)


def validate(predictions: Sequence[Prediction]) -> ValidationReport:
    """Check a submission against every rule we know of.

    Errors block submission; warnings are judgement calls worth a second look.
    """
    report = ValidationReport()

    if not predictions:
        report.errors.append("submission is empty")
        return report

    if len(predictions) > MAX_ROWS:
        report.errors.append(f"{len(predictions)} rows exceeds the {MAX_ROWS}-row limit")

    proband_ids = {p.proband_id for p in predictions}
    if len(proband_ids) > 1:
        report.errors.append(f"multiple proband_ids in one file: {sorted(proband_ids)}")
    if any(not p.proband_id for p in predictions):
        report.errors.append("proband_id is blank in at least one row")
    wrong = sorted(pid for pid in proband_ids if pid and pid != PROBAND_ID)
    if wrong:
        report.errors.append(
            f"proband_id must be {PROBAND_ID!r}, got {wrong}. The Space rejects anything else "
            "with 'Unknown proband_id'; the VCF sample name is not the proband id."
        )

    # --- EPCR discipline (see the organizers' published rules) ---------
    epcrs = [float(p.epcr) for p in predictions]
    duplicates = {value for value in epcrs if epcrs.count(value) > 1}
    if duplicates:
        report.errors.append(
            f"duplicate epcr values {sorted(duplicates)} — ties cannot be separated by any "
            "threshold and permanently depress F-max; make every epcr distinct"
        )

    primaries = [p for p in predictions if p.finding_type == "primary"]
    secondaries = [p for p in predictions if p.finding_type == "secondary"]
    if not primaries:
        report.errors.append("no row is marked finding_type='primary'")
    if secondaries and primaries:
        highest_secondary = max(float(p.epcr) for p in secondaries)
        lowest_primary = min(float(p.epcr) for p in primaries)
        if highest_secondary >= lowest_primary:
            report.warnings.append(
                f"a secondary finding (epcr={highest_secondary:g}) ranks at or above a primary "
                f"(epcr={lowest_primary:g}); secondary findings should sit in the lowest band. "
                "The FAQ says they cannot hurt your score, but the organizers' scorer never "
                "reads finding_type - a secondary row above a primary takes its rank and adds "
                "false positives to the F-max sweep"
            )

    top = max(predictions, key=lambda p: float(p.epcr))
    if float(top.epcr) < 0.5:
        report.warnings.append(
            f"top epcr is only {float(top.epcr):g}; if you believe a candidate, say so — "
            "relative ranking is what is scored, but a flat distribution wastes F-max headroom"
        )

    # --- Compound-het representation ----------------------------------------
    # The risk: a pair row and a standalone row of one of its alleles sitting at
    # indistinguishable confidence, so the F-max sweep cannot separate them and
    # the pair reads as two single guesses.
    #
    # 🔴 This used to test ``abs(a - b) < 0.1``, an ABSOLUTE epsilon, which is
    # scale-dependent and therefore backwards. Once ``epcr`` values are well
    # calibrated, the hedges are small (0.06, 0.04, 0.02) and every neighbouring
    # pair falls inside 0.1 — so the rule fired on a CORRECT submission. Inflate
    # the same ladder to 0.25/0.13 and it goes quiet. A check that rewards
    # overconfidence is worse than no check.
    #
    # The scale-free question is whether a threshold sweep can tell the two rows
    # apart, so compare RELATIVE separation.
    singles: dict[Variant, float] = {p.variant_1: float(p.epcr) for p in predictions if not p.is_pair}
    for pair in (p for p in predictions if p.is_pair):
        assert pair.variant_2 is not None
        for allele in (pair.variant_1, pair.variant_2):
            if allele not in singles:
                continue
            single, paired = singles[allele], float(pair.epcr)
            scale = max(single, paired)
            if scale > 0 and abs(single - paired) / scale < PAIR_SINGLE_RELATIVE_EPSILON:
                report.warnings.append(
                    f"{allele} appears both inside a compound-het pair (epcr {paired:g}) and as a "
                    f"standalone row (epcr {single:g}) at within {PAIR_SINGLE_RELATIVE_EPSILON:.0%} of the same "
                    "confidence; separate them or the pair may be scored as two single guesses"
                )

    # --- Duplicate hypotheses ------------------------------------------------
    seen: set[tuple[Variant, Variant | None]] = set()
    for prediction in predictions:
        key = (prediction.variant_1, prediction.variant_2)
        reverse = (prediction.variant_2, prediction.variant_1) if prediction.is_pair else key
        if key in seen or reverse in seen:
            report.errors.append(f"duplicate prediction: {prediction.variant_1} / {prediction.variant_2}")
        seen.add(key)

    return report


def write(predictions: Sequence[Prediction], path: str | Path, *, force: bool = False) -> Path:
    """Validate and write a Track 1 predictions CSV.

    Rows are written in descending ``epcr`` order, so rank order is visible in
    the file itself.

    Args:
        predictions: The rows to write.
        path: Destination CSV path.
        force: Write even if validation produced errors. For testing only.

    Raises:
        SubmissionError: if validation fails and ``force`` is not set.
    """
    report = validate(predictions)
    if not report.ok and not force:
        raise SubmissionError("submission failed validation:\n" + report.render())

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(predictions, key=lambda p: float(p.epcr), reverse=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(COLUMNS))
        writer.writeheader()
        for prediction in ordered:
            writer.writerow(prediction.to_row())
    return destination


def read(path: str | Path) -> list[Prediction]:
    """Read a predictions CSV back into :class:`Prediction` objects.

    Accepts either contig convention on input; the objects normalise internally.
    """
    predictions: list[Prediction] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = set(COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise SubmissionError(f"CSV is missing required columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                variant_2 = None
                if (row.get("chrom_2") or "").strip():
                    variant_2 = Variant(
                        contig=row["chrom_2"],
                        position=int(row["pos_2"]),
                        ref=row["ref_2"],
                        alt=row["alt_2"],
                    )
                predictions.append(
                    Prediction(
                        proband_id=(row.get("proband_id") or "").strip(),
                        variant_1=Variant(
                            contig=row["chrom_1"],
                            position=int(row["pos_1"]),
                            ref=row["ref_1"],
                            alt=row["alt_1"],
                        ),
                        variant_2=variant_2,
                        epcr=float(row["epcr"]),
                        finding_type=(row.get("finding_type") or "primary").strip().lower(),
                        notes=(row.get("notes") or "").strip(),
                    )
                )
            except (CoordinateError, SubmissionError, ValueError, KeyError) as exc:
                raise SubmissionError(f"line {line_number}: {exc}") from exc
    return predictions


def check_file(path: str | Path) -> ValidationReport:
    """Read a CSV and validate it, reporting parse failures as errors."""
    try:
        predictions = read(path)
    except SubmissionError as exc:
        return ValidationReport(errors=[str(exc)])
    report = validate(predictions)

    # A file-level check that only applies to on-disk files: the contig column
    # must already carry the chr prefix, because that is what gets uploaded.
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.DictReader(handle), start=2):
            for column in ("chrom_1", "chrom_2"):
                value = (row.get(column) or "").strip()
                if value and not value.lower().startswith("chr"):
                    report.errors.append(
                        f"line {line_number}: {column}={value!r} is missing the 'chr' prefix required "
                        "by the submission template (the source VCF uses bare contigs — convert on write)"
                    )
    return report


def suggested_filename(username: str, approach: str, *, extension: str = "csv") -> str:
    """Build a filename in the organizers' recommended form.

    Examples:
        >>> suggested_filename("jane-doe", "Blind Exomiser + ACMG")
        'jane-doe_blind-exomiser-acmg.csv'
    """
    import re

    slug = re.sub(r"[^a-z0-9]+", "-", approach.lower()).strip("-")
    return f"{username}_{slug}.{extension}"


def format_table(predictions: Iterable[Prediction]) -> str:
    """Render predictions as a readable table for review before submission."""
    rows = sorted(predictions, key=lambda p: float(p.epcr), reverse=True)
    lines = [f"{'rank':>4}  {'epcr':>6}  {'type':<9}  variant(s)"]
    lines.append("-" * 72)
    for rank, prediction in enumerate(rows, start=1):
        variants = str(prediction.variant_1)
        if prediction.variant_2 is not None:
            variants += f"  +  {prediction.variant_2}"
        lines.append(f"{rank:>4}  {float(prediction.epcr):>6.3f}  {prediction.finding_type:<9}  {variants}")
        if prediction.notes:
            lines.append(f"{'':>4}  {'':>6}  {'':<9}  ↳ {prediction.notes}")
    return "\n".join(lines)
