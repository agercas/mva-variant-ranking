"""Contig-naming and variant-representation helpers.

🔴 **The single highest-risk detail in Track 1.**

The hackathon VCF is aligned to
``GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta``
and therefore uses **bare contig names** (``1``, ``2``, … ``X``, ``Y``, ``MT``).
The Track 1 submission template requires the **``chr`` prefix** (``chr15``).

Two leaderboard entries scored 0.0 — exactly what a convention mismatch
produces. Everything crossing that boundary goes through this module.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

#: Valid contigs, in the order a human expects to see them.
AUTOSOMES: tuple[str, ...] = tuple(str(i) for i in range(1, 23))
ALLOSOMES: tuple[str, ...] = ("X", "Y")
MITO: tuple[str, ...] = ("M", "MT")
CANONICAL_CONTIGS: frozenset[str] = frozenset(AUTOSOMES + ALLOSOMES + MITO)

_BASES = re.compile(r"^[ACGTN]+$", re.IGNORECASE)


class CoordinateError(ValueError):
    """Raised when a contig, position or allele cannot be interpreted."""


def normalise_contig(contig: str) -> str:
    """Strip any ``chr`` prefix and canonicalise case.

    Args:
        contig: A contig name in any common convention (``chr15``, ``15``,
            ``CHR15``, ``chrX``, ``chrM``).

    Returns:
        The bare contig name, upper-cased for sex and mitochondrial contigs.

    Raises:
        CoordinateError: if the contig is not a canonical human contig.

    Examples:
        >>> normalise_contig("chr15")
        '15'
        >>> normalise_contig("X")
        'X'
    """
    bare = str(contig).strip()
    if bare.lower().startswith("chr"):
        bare = bare[3:]
    if bare.upper() in {*ALLOSOMES, *MITO}:
        bare = bare.upper()
    if bare not in CANONICAL_CONTIGS:
        raise CoordinateError(
            f"{contig!r} is not a canonical human contig; "
            f"expected one of 1-22, X, Y, M/MT (with or without a 'chr' prefix)"
        )
    return bare


def to_vcf_contig(contig: str) -> str:
    """Return the bare form used by this dataset's VCF (``15``).

    Examples:
        >>> to_vcf_contig("chr15")
        '15'
    """
    return normalise_contig(contig)


def to_submission_contig(contig: str) -> str:
    """Return the ``chr``-prefixed form required by the submission CSV (``chr15``).

    Examples:
        >>> to_submission_contig("15")
        'chr15'
        >>> to_submission_contig("chrX")
        'chrX'
    """
    return f"chr{normalise_contig(contig)}"


def validate_allele(allele: str, *, field: str = "allele") -> str:
    """Check that an allele is a non-empty run of ACGTN and return it upper-cased.

    Raises:
        CoordinateError: if the allele is empty or contains non-nucleotide
            characters. Symbolic alleles (``<DEL>``) and breakend notation are
            rejected: the submission schema has no representation for them.
    """
    value = str(allele).strip().upper()
    if not value:
        raise CoordinateError(f"{field} is empty")
    if not _BASES.match(value):
        raise CoordinateError(
            f"{field} {allele!r} is not a plain nucleotide sequence; "
            "symbolic and breakend alleles cannot be represented in the Track 1 schema"
        )
    return value


@dataclass(frozen=True, order=True)
class Variant:
    """A single genomic variant, stored in the dataset's bare-contig convention.

    Construct with any contig convention; the stored form is always bare, so
    equality and hashing are convention-independent — ``chr15:40163000C>T`` and
    ``15:40163000C>T`` are the same :class:`Variant`.
    """

    contig: str
    position: int
    ref: str
    alt: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "contig", normalise_contig(self.contig))
        object.__setattr__(self, "ref", validate_allele(self.ref, field="ref"))
        object.__setattr__(self, "alt", validate_allele(self.alt, field="alt"))
        if int(self.position) < 1:
            raise CoordinateError(f"position must be 1-based and positive, got {self.position}")
        object.__setattr__(self, "position", int(self.position))
        if self.ref == self.alt:
            raise CoordinateError(f"ref and alt are identical ({self.ref}) — not a variant")

    @classmethod
    def parse(cls, text: str) -> Variant:
        """Parse ``"chr15-40163000-C-T"`` or ``"15:40163000:C:T"``.

        Examples:
            >>> Variant.parse("chr15-40163000-C-T")
            Variant(contig='15', position=40163000, ref='C', alt='T')
        """
        parts = re.split(r"[-:]", text.strip())
        if len(parts) != 4:
            raise CoordinateError(
                f"cannot parse variant {text!r}; expected CONTIG-POS-REF-ALT or CONTIG:POS:REF:ALT"
            )
        contig, position, ref, alt = parts
        return cls(contig=contig, position=int(position), ref=ref, alt=alt)

    @property
    def submission_contig(self) -> str:
        """The ``chr``-prefixed contig required by the submission CSV."""
        return to_submission_contig(self.contig)

    @property
    def is_snv(self) -> bool:
        """True for a single-nucleotide substitution."""
        return len(self.ref) == 1 and len(self.alt) == 1

    def to_gnomad_id(self) -> str:
        """gnomAD's variant id form, ``"15-40163000-C-T"`` (no ``chr``)."""
        return f"{self.contig}-{self.position}-{self.ref}-{self.alt}"

    def to_myvariant_id(self) -> str:
        """MyVariant.info's hgvs-style id, ``"chr15:g.40163000C>T"`` (with ``chr``)."""
        if self.is_snv:
            return f"chr{self.contig}:g.{self.position}{self.ref}>{self.alt}"
        # Indels use delins notation spanning the reference allele.
        end = self.position + len(self.ref) - 1
        return f"chr{self.contig}:g.{self.position}_{end}delins{self.alt}"

    def to_spdi(self) -> str:
        """NCBI SPDI-ish form used by several APIs."""
        return f"{self.contig}:{self.position - 1}:{self.ref}:{self.alt}"

    def __str__(self) -> str:
        return f"{self.contig}:{self.position}{self.ref}>{self.alt}"


def to_gnomad_id(contig: str, position: int, ref: str, alt: str) -> str:
    """Convenience wrapper producing a gnomAD variant id from loose parts."""
    return Variant(contig=contig, position=position, ref=ref, alt=alt).to_gnomad_id()
