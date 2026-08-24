"""Cohort, donor, preprocessing, and locus leakage controls."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
from math import isfinite
from typing import Any, Hashable, Iterable, Mapping, Sequence


class LeakageError(RuntimeError):
    """Raised when an identity or preprocessing object crosses a split."""


CHAMPION_ELIGIBLE_EXPOSURES = frozenset({"clean_declared", "target_label_unexposed"})


def salted_signature(value: bytes | str, *, salt: bytes) -> str:
    """Create a non-reversible stored signature; raw donor data are not retained."""

    if len(salt) < 16:
        raise LeakageError("fingerprint salt must contain at least 128 bits")
    payload = value.encode("utf-8") if isinstance(value, str) else value
    return hashlib.sha256(salt + b"\0" + payload).hexdigest()


def assert_checkpoint_champion_eligible(model: Mapping[str, Any]) -> None:
    exposure = str(model.get("exposure_status", "unknown"))
    license_status = str(model.get("license_status", "unknown"))
    if exposure not in CHAMPION_ELIGIBLE_EXPOSURES:
        raise LeakageError(f"checkpoint exposure is not champion-eligible: {exposure}")
    if license_status not in {
        "open",
        "open_source",
        "permissive",
        "redistributable",
        "project_owned",
        "MIT",
        "BSD-3-Clause",
        "Apache-2.0",
        "CC-BY-4.0",
    }:
        raise LeakageError(f"checkpoint license is not open-champion eligible: {license_status}")
    if not model.get("checkpoint_revision") or not model.get("checkpoint_sha256"):
        raise LeakageError("checkpoint revision and SHA-256 must be exact")


def assert_group_integrity(
    rows: Iterable[Mapping[str, Any]],
    *,
    group_keys: Sequence[str],
    fold_key: str = "outer_fold",
) -> None:
    """Require every donor/cohort alias group to occupy exactly one outer fold."""

    memberships: dict[tuple[Hashable, ...], set[Hashable]] = defaultdict(set)
    for row in rows:
        try:
            group = tuple(row[key] for key in group_keys)
            fold = row[fold_key]
        except KeyError as error:
            raise LeakageError(f"split row is missing required key: {error.args[0]}") from error
        memberships[group].add(fold)
    leaking = [group for group, folds in memberships.items() if len(folds) > 1]
    if leaking:
        raise LeakageError(f"{len(leaking)} identity groups cross outer folds")


def assert_preprocessor_fit_scope(
    provenance: Mapping[str, Any],
    *,
    allowed_training_ids: set[str],
    sealed_ids: set[str],
) -> None:
    fit_ids = {str(item) for item in provenance.get("fit_partition_ids", [])}
    if not fit_ids:
        raise LeakageError("preprocessor provenance has no fit_partition_ids")
    if fit_ids.intersection(sealed_ids):
        raise LeakageError("preprocessing was fit using sealed observations")
    unexpected = fit_ids.difference(allowed_training_ids)
    if unexpected:
        raise LeakageError(
            "preprocessing was fit outside the training partition: " + ", ".join(sorted(unexpected))
        )


@dataclass(frozen=True, slots=True)
class Variant:
    variant_id: str
    chromosome: str
    position: int
    gene_ids: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class SplitWindow:
    window_id: str
    chromosome: str
    start: int
    end: int
    outer_role: str


def assert_sequence_window_isolation(
    windows: Sequence[SplitWindow], *, boundary_buffer_bp: int
) -> None:
    """Reject cross-role windows separated by less than the frozen receptive-field buffer."""

    if (
        isinstance(boundary_buffer_bp, bool)
        or not isinstance(boundary_buffer_bp, int)
        or boundary_buffer_bp < 0
    ):
        raise LeakageError("boundary_buffer_bp must be a nonnegative integer")
    identifiers = [window.window_id for window in windows]
    if len(identifiers) != len(set(identifiers)):
        raise LeakageError("sequence window identifiers must be unique")
    by_chromosome: dict[str, list[SplitWindow]] = defaultdict(list)
    for window in windows:
        if not window.chromosome or not window.outer_role:
            raise LeakageError("sequence windows require chromosome and outer_role")
        if (
            isinstance(window.start, bool)
            or isinstance(window.end, bool)
            or not isinstance(window.start, int)
            or not isinstance(window.end, int)
            or window.start < 0
            or window.end <= window.start
        ):
            raise LeakageError(f"invalid sequence window coordinates: {window.window_id}")
        by_chromosome[window.chromosome].append(window)
    leaking: list[tuple[str, str]] = []
    for chromosome_windows in by_chromosome.values():
        ordered = sorted(
            chromosome_windows, key=lambda item: (item.start, item.end, item.window_id)
        )
        active: list[SplitWindow] = []
        for current in ordered:
            active = [
                previous
                for previous in active
                if previous.end + boundary_buffer_bp > current.start
            ]
            for previous in active:
                if previous.outer_role != current.outer_role:
                    leaking.append((previous.window_id, current.window_id))
                    if len(leaking) == 10:
                        break
            if len(leaking) == 10:
                break
            active.append(current)
        if len(leaking) == 10:
            break
    if leaking:
        pairs = ", ".join(f"{left}/{right}" for left, right in leaking)
        raise LeakageError(
            "sequence windows cross outer roles within the receptive-field buffer: "
            + pairs
        )


class _UnionFind:
    def __init__(self, values: Iterable[str]):
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[max(left_root, right_root)] = min(left_root, right_root)


def merged_variant_blocks(
    variants: Sequence[Variant],
    *,
    ld_edges_by_ancestry: Mapping[str, Iterable[tuple[str, str, float]]],
    distance_bp: int = 1_000_000,
    ld_r2_threshold: float = 0.8,
) -> dict[str, str]:
    """Merge variants connected by distance, gene, or LD in any ancestry."""

    if (
        isinstance(distance_bp, bool)
        or not isinstance(distance_bp, int)
        or distance_bp < 0
    ):
        raise LeakageError("distance_bp must be a nonnegative integer")
    if isinstance(ld_r2_threshold, bool):
        raise LeakageError("ld_r2_threshold must be numeric")
    try:
        threshold = float(ld_r2_threshold)
    except (TypeError, ValueError) as error:
        raise LeakageError("ld_r2_threshold must be numeric") from error
    if not isfinite(threshold) or not 0.0 <= threshold <= 1.0:
        raise LeakageError("ld_r2_threshold must be finite and between zero and one")
    if any(not isinstance(variant, Variant) for variant in variants):
        raise LeakageError("variants must contain Variant records")
    identifiers = [variant.variant_id for variant in variants]
    if len(identifiers) != len(set(identifiers)):
        raise LeakageError("variant identifiers must be unique")
    for variant in variants:
        if (
            not isinstance(variant.variant_id, str)
            or not variant.variant_id
            or not isinstance(variant.chromosome, str)
            or not variant.chromosome
        ):
            raise LeakageError("variants require nonempty identifiers and chromosomes")
        if (
            isinstance(variant.position, bool)
            or not isinstance(variant.position, int)
            or variant.position < 1
        ):
            raise LeakageError("variant positions must be positive one-based integers")
        if (
            any(not isinstance(gene_id, str) or not gene_id for gene_id in variant.gene_ids)
            or len(set(variant.gene_ids)) != len(variant.gene_ids)
        ):
            raise LeakageError("variant gene_ids must be unique nonempty identifiers")
    union = _UnionFind(identifiers)
    by_chromosome: dict[str, list[Variant]] = defaultdict(list)
    by_gene: dict[str, list[str]] = defaultdict(list)
    for variant in variants:
        by_chromosome[variant.chromosome].append(variant)
        for gene_id in variant.gene_ids:
            by_gene[gene_id].append(variant.variant_id)
    for chromosome_variants in by_chromosome.values():
        ordered = sorted(chromosome_variants, key=lambda item: item.position)
        left = 0
        for right, variant in enumerate(ordered):
            while ordered[left].position < variant.position - distance_bp:
                left += 1
            for neighbor in ordered[left:right]:
                union.union(variant.variant_id, neighbor.variant_id)
    for members in by_gene.values():
        for member in members[1:]:
            union.union(members[0], member)
    known = set(identifiers)
    for ancestry, edges in ld_edges_by_ancestry.items():
        if not isinstance(ancestry, str) or not ancestry:
            raise LeakageError("LD ancestry identifiers must be nonempty strings")
        for edge in edges:
            if not isinstance(edge, (tuple, list)) or len(edge) != 3:
                raise LeakageError(f"{ancestry} LD edges must contain left, right, and r2")
            left, right, r_squared = edge
            if left not in known or right not in known:
                raise LeakageError(f"{ancestry} LD edge names an unknown variant")
            if isinstance(r_squared, bool):
                raise LeakageError(f"{ancestry} LD r2 must be numeric")
            try:
                checked_r_squared = float(r_squared)
            except (TypeError, ValueError) as error:
                raise LeakageError(f"{ancestry} LD r2 must be numeric") from error
            if not isfinite(checked_r_squared) or not 0.0 <= checked_r_squared <= 1.0:
                raise LeakageError(f"{ancestry} LD r2 must be finite and within [0, 1]")
            if checked_r_squared >= threshold:
                union.union(left, right)
    roots = {identifier: union.find(identifier) for identifier in identifiers}
    stable_roots = {root: index for index, root in enumerate(sorted(set(roots.values())), start=1)}
    return {
        identifier: f"variant-block-{stable_roots[root]:06d}"
        for identifier, root in roots.items()
    }
