"""Reference-safe coordinate and allele primitives for sequence adapters."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from .reference import CoordinateContract, ReferenceError


_DNA = frozenset("ACGTN")
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _dna(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ReferenceError(f"{field_name} must be a non-empty DNA string")
    normalized = value.upper()
    invalid = sorted(set(normalized) - _DNA)
    if invalid:
        raise ReferenceError(
            f"{field_name} contains unsupported base(s): {', '.join(invalid)}"
        )
    return normalized


def reverse_complement(sequence: str) -> str:
    """Return an uppercase reverse complement over the explicit A/C/G/T/N alphabet."""

    return _dna(sequence, "sequence").translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True, slots=True)
class BedInterval:
    contig: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


def validate_bed_interval(
    contig: str,
    start: int,
    end: int,
    *,
    coordinate_contract: CoordinateContract,
    contig_lengths: Mapping[str, int] | None = None,
) -> BedInterval:
    """Validate one 0-based, half-open interval against the admitted contig set."""

    normalized_contig = coordinate_contract.require_analysis_contig(contig)
    if isinstance(start, bool) or not isinstance(start, int):
        raise ReferenceError("BED start must be an integer")
    if isinstance(end, bool) or not isinstance(end, int):
        raise ReferenceError("BED end must be an integer")
    if start < 0 or end <= start:
        raise ReferenceError("BED interval must satisfy 0 <= start < end")
    if contig_lengths is not None:
        length = contig_lengths.get(normalized_contig)
        if isinstance(length, bool) or not isinstance(length, int) or length <= 0:
            raise ReferenceError(
                f"positive contig length is required for {normalized_contig}"
            )
        if end > length:
            raise ReferenceError(
                f"BED interval crosses the {normalized_contig} boundary: {end} > {length}"
            )
    return BedInterval(normalized_contig, start, end)


def vcf_to_bed(
    contig: str,
    position_1based: int,
    ref: str,
    alt: str,
    *,
    coordinate_contract: CoordinateContract,
    contig_lengths: Mapping[str, int] | None = None,
) -> BedInterval:
    """Convert a simple VCF allele record to its reference-spanning BED interval."""

    if isinstance(position_1based, bool) or not isinstance(position_1based, int):
        raise ReferenceError("VCF position must be an integer")
    if position_1based < 1:
        raise ReferenceError("VCF position must be one-based and positive")
    normalized_ref = _dna(ref, "REF")
    _dna(alt, "ALT")
    start = position_1based - 1
    return validate_bed_interval(
        contig,
        start,
        start + len(normalized_ref),
        coordinate_contract=coordinate_contract,
        contig_lengths=contig_lengths,
    )


def require_ref_allele(observed_reference: str, declared_ref: str) -> str:
    """Reject a strand/build/coordinate mismatch before variant scoring."""

    observed = _dna(observed_reference, "observed_reference")
    declared = _dna(declared_ref, "REF")
    if observed != declared:
        raise ReferenceError(
            f"REF allele mismatch: declared {declared}, reference contains {observed}"
        )
    return declared


def centered_window(
    contig: str,
    anchor_0based: int,
    width: int,
    *,
    coordinate_contract: CoordinateContract,
    contig_lengths: Mapping[str, int],
) -> BedInterval:
    """Create an exact-width half-open window anchored on a zero-based variant base."""

    if isinstance(anchor_0based, bool) or not isinstance(anchor_0based, int):
        raise ReferenceError("window anchor must be an integer")
    if isinstance(width, bool) or not isinstance(width, int) or width <= 0:
        raise ReferenceError("window width must be a positive integer")
    start = anchor_0based - width // 2
    return validate_bed_interval(
        contig,
        start,
        start + width,
        coordinate_contract=coordinate_contract,
        contig_lengths=contig_lengths,
    )


__all__ = [
    "BedInterval",
    "centered_window",
    "require_ref_allele",
    "reverse_complement",
    "validate_bed_interval",
    "vcf_to_bed",
]
