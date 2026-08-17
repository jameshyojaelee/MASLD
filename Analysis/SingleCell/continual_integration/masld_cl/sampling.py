"""Seeded donor-balanced optimization, replay, Fisher, and BI sampling."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from typing import Hashable

import numpy as np
import torch


class SamplingError(ValueError):
    pass


def capped_group_indices(
    groups: Sequence[tuple[Hashable, ...]], cap: int, seed: int
) -> np.ndarray:
    if cap <= 0:
        raise SamplingError("cap must be positive")
    members: dict[tuple[Hashable, ...], list[int]] = defaultdict(list)
    for index, group in enumerate(groups):
        members[tuple(group)].append(index)
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    for group in sorted(members, key=lambda x: tuple(map(str, x))):
        values = np.asarray(members[group], dtype=np.int64)
        if len(values) > cap:
            values = rng.choice(values, size=cap, replace=False)
        selected.extend(int(x) for x in values)
    return np.asarray(sorted(selected), dtype=np.int64)


def stratified_fraction_indices(
    candidate_indices: Sequence[int],
    strata: Sequence[tuple[Hashable, ...]],
    fraction: float,
    seed: int,
) -> np.ndarray:
    """Sample the requested fraction within every nonempty stratum.

    Rounding uses nearest integer, with one item retained for nonzero fractions.
    The same source index cannot be emitted twice.
    """
    if not 0 <= fraction <= 1:
        raise SamplingError("fraction must be between zero and one")
    if len(candidate_indices) != len(strata):
        raise SamplingError("candidate_indices and strata must have equal length")
    by_group: dict[tuple[Hashable, ...], list[int]] = defaultdict(list)
    for index, group in zip(candidate_indices, strata):
        by_group[tuple(group)].append(int(index))
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    for group in sorted(by_group, key=lambda x: tuple(map(str, x))):
        values = np.asarray(sorted(set(by_group[group])), dtype=np.int64)
        if len(values) != len(by_group[group]):
            raise SamplingError(f"duplicate source index in stratum {group}")
        target = 0 if fraction == 0 else max(1, int(np.floor(len(values) * fraction + 0.5)))
        target = min(target, len(values))
        if target:
            selected.extend(int(x) for x in rng.choice(values, target, replace=False))
    result = np.asarray(sorted(selected), dtype=np.int64)
    if len(result) != len(np.unique(result)):
        raise SamplingError("sampler emitted a duplicate index")
    return result


def balanced_exact_fraction_indices(
    candidate_indices: Sequence[int],
    strata: Sequence[tuple[Hashable, ...]],
    fraction: float,
    seed: int,
) -> np.ndarray:
    """Sample an exact fraction with max-min allocation across strata.

    Small strata saturate first. Every unsaturated donor-lineage-preparation
    stratum then differs by at most one selected cell. The total is exactly the
    nearest integer fraction of all eligible candidates.
    """
    candidates = np.asarray(candidate_indices, dtype=np.int64)
    if len(candidates) != len(strata) or len(np.unique(candidates)) != len(candidates):
        raise SamplingError("balanced candidates must be unique and match strata")
    if not 0 <= fraction <= 1:
        raise SamplingError("fraction must be between zero and one")
    target = int(np.floor(len(candidates) * fraction + 0.5))
    if target == 0:
        return np.asarray([], dtype=np.int64)
    groups: dict[tuple[Hashable, ...], list[int]] = defaultdict(list)
    for local, stratum in enumerate(strata):
        groups[tuple(stratum)].append(local)
    ordered = sorted(groups, key=lambda value: tuple(map(str, value)))
    capacities = {group: len(groups[group]) for group in ordered}
    allocation = {group: 0 for group in ordered}
    remaining = target
    level = 0
    while remaining:
        active = [group for group in ordered if capacities[group] > level]
        if not active:
            raise SamplingError("balanced allocation exhausted before reaching target")
        if remaining >= len(active):
            for group in active:
                allocation[group] += 1
            remaining -= len(active)
            level += 1
        else:
            # Seeded rotation avoids always awarding the remainder to lexical-first strata.
            rng = np.random.default_rng(seed + 7919 * (level + 1))
            chosen = rng.choice(len(active), size=remaining, replace=False)
            for index in sorted(map(int, chosen)):
                allocation[active[index]] += 1
            remaining = 0
    rng = np.random.default_rng(seed)
    selected = []
    for group in ordered:
        local = np.asarray(groups[group], dtype=np.int64)
        count = allocation[group]
        if count:
            selected.extend(candidates[rng.choice(local, count, replace=False)])
    result = np.asarray(sorted(map(int, selected)), dtype=np.int64)
    if len(result) != target or len(np.unique(result)) != target:
        raise SamplingError("balanced sampler did not return the exact unique target")
    unsaturated = [allocation[group] for group in ordered if allocation[group] < capacities[group]]
    if unsaturated and max(unsaturated) - min(unsaturated) > 1:
        raise SamplingError("balanced sampler violated max-min allocation")
    return result


def concatenate_query_and_replay(query_indices: Sequence[int], replay_indices: Sequence[int]) -> np.ndarray:
    query = np.asarray(query_indices, dtype=np.int64)
    replay = np.asarray(replay_indices, dtype=np.int64)
    if len(np.unique(query)) != len(query) or len(np.unique(replay)) != len(replay):
        raise SamplingError("query and replay inputs must each contain unique cells")
    overlap = np.intersect1d(query, replay)
    if len(overlap):
        raise SamplingError("query and replay cells overlap")
    result = np.concatenate([query, replay])
    if len(np.unique(result)) != len(result):
        raise SamplingError("training object contains a cell more than once")
    return result


def donor_train_validation_split(
    donor_ids: Sequence[str], validation_fraction: float, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    if not 0 < validation_fraction < 1:
        raise SamplingError("validation_fraction must be strictly between zero and one")
    donors = np.asarray(donor_ids, dtype=object)
    unique = np.asarray(sorted(set(donors)), dtype=object)
    if len(unique) < 2:
        raise SamplingError("at least two biological donors are required")
    rng = np.random.default_rng(seed)
    n_validation = max(1, int(np.floor(len(unique) * validation_fraction + 0.5)))
    n_validation = min(n_validation, len(unique) - 1)
    validation_donors = set(rng.choice(unique, n_validation, replace=False))
    validation = np.flatnonzero(np.fromiter((x in validation_donors for x in donors), bool))
    train = np.flatnonzero(np.fromiter((x not in validation_donors for x in donors), bool))
    if set(donors[train]) & set(donors[validation]):
        raise SamplingError("donor leakage between train and validation")
    return train.astype(np.int64), validation.astype(np.int64)


def masked_augmentations(
    expression: torch.Tensor,
    *,
    n_augmentations: int = 200,
    mask_fraction: float = 0.50,
    generator: torch.Generator | None = None,
) -> torch.Tensor:
    """Return TTA inputs with exactly the requested fraction of genes zeroed."""
    if expression.ndim != 2:
        raise SamplingError("expression must have shape cells x genes")
    if not 0 < mask_fraction < 1 or n_augmentations <= 0:
        raise SamplingError("invalid augmentation parameters")
    n_genes = expression.shape[1]
    n_masked = int(round(n_genes * mask_fraction))
    outputs = []
    for _ in range(n_augmentations):
        order = torch.randperm(n_genes, generator=generator, device=expression.device)
        mask = torch.ones(n_genes, dtype=expression.dtype, device=expression.device)
        mask[order[:n_masked]] = 0
        outputs.append(expression * mask)
    return torch.stack(outputs)


def bregman_information_lse(latents: torch.Tensor) -> torch.Tensor:
    """Bregman information for LogSumExp over TTA latent representations."""
    if latents.ndim != 3:
        raise SamplingError("latents must have shape augmentations x cells x latent")
    expected_lse = torch.logsumexp(latents, dim=-1).mean(dim=0)
    lse_expected = torch.logsumexp(latents.mean(dim=0), dim=-1)
    return expected_lse - lse_expected


def exact_quantile_replay(scores: Sequence[float], size: int, mode: str) -> np.ndarray:
    """Select exactly `size` cells by deterministic BI ranking."""
    values = np.asarray(scores, dtype=float)
    if not np.all(np.isfinite(values)) or not 0 <= size <= len(values):
        raise SamplingError("invalid BI scores or replay size")
    stable = np.lexsort((np.arange(len(values)), values))
    if mode == "bottom":
        result = stable[:size]
    elif mode == "top":
        result = stable[::-1][:size]
    elif mode == "step":
        # Evenly cover the full ordered distribution, including both endpoints.
        if size == 0:
            result = np.asarray([], dtype=np.int64)
        elif size == 1:
            result = stable[[len(stable) // 2]]
        else:
            positions = np.rint(np.linspace(0, len(stable) - 1, size)).astype(int)
            result = stable[positions]
    else:
        raise SamplingError("mode must be bottom, top, or step")
    if len(result) != size or len(np.unique(result)) != size:
        raise SamplingError("BI selector did not return the exact unique buffer size")
    return np.asarray(result, dtype=np.int64)


def stratified_exact_bi_replay(
    candidate_indices: Sequence[int], scores: Sequence[float],
    strata: Sequence[tuple[Hashable, ...]], size: int, mode: str,
    target_counts: Mapping[tuple[Hashable, ...], int] | None = None,
) -> np.ndarray:
    """Allocate an exact BI buffer proportionally, then rank within each stratum."""
    candidates = np.asarray(candidate_indices, dtype=np.int64)
    scores = np.asarray(scores, dtype=float)
    if len(candidates) != len(scores) or len(candidates) != len(strata):
        raise SamplingError("BI candidates, scores, and strata must have equal length")
    if len(np.unique(candidates)) != len(candidates) or not 0 <= size <= len(candidates):
        raise SamplingError("BI candidates must be unique and buffer size valid")
    groups: dict[tuple[Hashable, ...], list[int]] = defaultdict(list)
    for local, stratum in enumerate(strata):
        groups[tuple(stratum)].append(local)
    ordered = sorted(groups, key=lambda x: tuple(map(str, x)))
    if target_counts is None:
        ideals = {
            group: size * len(groups[group]) / len(candidates) for group in ordered
        }
        allocation = {group: int(np.floor(ideals[group])) for group in ordered}
        remaining = size - sum(allocation.values())
        priority = sorted(
            ordered,
            key=lambda group: (
                -(ideals[group] - allocation[group]), tuple(map(str, group))
            ),
        )
        for group in priority[:remaining]:
            allocation[group] += 1
    else:
        normalized = {tuple(group): int(count) for group, count in target_counts.items()}
        extra = set(normalized) - set(ordered)
        if extra:
            raise SamplingError(f"BI target counts contain unknown strata: {extra}")
        allocation = {group: normalized.get(group, 0) for group in ordered}
        if (
            any(count < 0 or count > len(groups[group]) for group, count in allocation.items())
            or sum(allocation.values()) != size
        ):
            raise SamplingError("BI target counts are invalid for the requested buffer")
    selected = []
    for group in ordered:
        local = np.asarray(groups[group], dtype=np.int64)
        ranked = exact_quantile_replay(scores[local], allocation[group], mode)
        selected.extend(candidates[local[ranked]])
    result = np.asarray(sorted(selected), dtype=np.int64)
    if len(result) != size or len(np.unique(result)) != size:
        raise SamplingError("stratified BI replay did not return the exact unique buffer size")
    return result
