#!/usr/bin/env python3
"""Pure genetics helpers used by the replay fixtures and lineage aggregator."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from atac_context_v3_lib import classify_genetic


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(value: str) -> str:
    return value.translate(COMPLEMENT)[::-1].upper()


@dataclass(frozen=True)
class NormalizedVariant:
    position_1based: int
    ref: str
    alt: str
    orientation: str


def _minimal_representation(position: int, ref: str, alt: str) -> tuple[int, str, str]:
    ref = ref.upper()
    alt = alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        position += 1
    return position, ref, alt


def normalize_against_reference(
    chrom: str,
    position_1based: int,
    allele1: str,
    allele2: str,
    fetch: Callable[[str, int, int], str],
) -> NormalizedVariant | None:
    """Orient, minimize, and left-align a biallelic variant on hg38.

    ``fetch`` uses zero-based half-open coordinates. One of the two supplied
    alleles must match the reference; no assumption is made about effect-allele
    ordering. Indels are explicitly left-aligned through their repeat context.
    """
    if not allele1 or not allele2:
        return None
    candidates = [
        (allele1.upper(), allele2.upper(), "forward"),
        (reverse_complement(allele1), reverse_complement(allele2), "reverse_complement"),
    ]
    for left, right, orientation in candidates:
        for ref, alt in ((left, right), (right, left)):
            observed = fetch(chrom, position_1based - 1, position_1based - 1 + len(ref)).upper()
            if observed != ref:
                continue
            pos, norm_ref, norm_alt = _minimal_representation(position_1based, ref, alt)
            if len(norm_ref) != len(norm_alt):
                while pos > 1:
                    previous = fetch(chrom, pos - 2, pos - 1).upper()
                    if not previous or norm_ref[-1] != previous or norm_alt[-1] != previous:
                        break
                    norm_ref = previous + norm_ref[:-1]
                    norm_alt = previous + norm_alt[:-1]
                    pos -= 1
                    pos, norm_ref, norm_alt = _minimal_representation(pos, norm_ref, norm_alt)
            final_ref = fetch(chrom, pos - 1, pos - 1 + len(norm_ref)).upper()
            if final_ref == norm_ref:
                return NormalizedVariant(pos, norm_ref, norm_alt, orientation)
    return None


def aggregate_lineage_mass(
    variants: list[dict[str, object]],
    overlap244: Callable[[str, int], bool],
    overlap281: Callable[[str, int], bool],
) -> dict[str, float | str]:
    mapped = sum(float(row["posterior"]) for row in variants if row.get("mapped"))
    mass244 = sum(
        float(row["posterior"])
        for row in variants
        if row.get("mapped") and overlap244(str(row["chrom"]), int(row["position_1based"]))
    )
    mass281 = sum(
        float(row["posterior"])
        for row in variants
        if row.get("mapped") and overlap281(str(row["chrom"]), int(row["position_1based"]))
    )
    shared = sum(
        float(row["posterior"])
        for row in variants
        if row.get("mapped")
        and overlap244(str(row["chrom"]), int(row["position_1based"]))
        and overlap281(str(row["chrom"]), int(row["position_1based"]))
    )
    return {
        "mapped_mass": mapped,
        "lost_mass": 1.0 - mapped,
        "gse244832_mass": mass244,
        "gse281367_mass": mass281,
        "shared_mass": shared,
        "state": classify_genetic(mapped, shared, mass244, mass281),
    }
