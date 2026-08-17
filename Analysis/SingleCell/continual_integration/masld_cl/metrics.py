"""Donor-level integration metrics and stratified bootstrap utilities."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np


class MetricError(ValueError):
    pass


def donor_centroids(latent: np.ndarray, donor_ids: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    latent = np.asarray(latent, dtype=float)
    donors = np.asarray(donor_ids, dtype=object)
    if latent.ndim != 2 or len(latent) != len(donors) or not np.isfinite(latent).all():
        raise MetricError("latent and donor IDs are invalid")
    unique = np.asarray(sorted(set(donors)), dtype=object)
    centroids = np.stack([latent[donors == donor].mean(axis=0) for donor in unique])
    return unique, centroids


def pairwise_euclidean(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    squared = (
        np.square(left).sum(axis=1, keepdims=True)
        + np.square(right).sum(axis=1)[None, :]
        - 2 * left @ right.T
    )
    return np.sqrt(np.maximum(squared, 0.0))


def mean_within_distance(values: np.ndarray) -> float:
    if len(values) < 2:
        raise MetricError("within-group distance requires at least two donors")
    distances = pairwise_euclidean(values, values)
    upper = distances[np.triu_indices(len(values), k=1)]
    return float(upper.mean())


def normalized_shift_control(reference: np.ndarray, controls: np.ndarray) -> float:
    """Group-centroid distance divided by mean within-group distance."""
    if len(reference) < 2 or len(controls) < 2:
        raise MetricError("shift-control requires at least two donors per group")
    cross = float(np.linalg.norm(reference.mean(axis=0) - controls.mean(axis=0)))
    denominator = 0.5 * (mean_within_distance(reference) + mean_within_distance(controls))
    if denominator <= 0:
        raise MetricError("shift-control denominator is zero")
    return cross / denominator


def standardized_case_control_separation(cases: np.ndarray, controls: np.ndarray) -> float:
    """Distance between donor-group means standardized by within-group scatter."""
    if len(cases) < 3 or len(controls) < 3:
        raise MetricError("disease separation requires at least three donors per group")
    case_center = cases.mean(axis=0)
    control_center = controls.mean(axis=0)
    between = float(np.linalg.norm(case_center - control_center))
    scatter = np.concatenate(
        [np.linalg.norm(cases - case_center, axis=1), np.linalg.norm(controls - control_center, axis=1)]
    )
    scale = float(np.sqrt(np.mean(np.square(scatter))))
    if scale <= 0:
        raise MetricError("within-group disease scale is zero")
    return between / scale


def _average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + end - 1)
        start = end
    return ranks


def spearman_correlation(left: Sequence[float], right: Sequence[float]) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape or left.ndim != 1 or len(left) < 3:
        raise MetricError("Spearman inputs must be equal vectors of length at least three")
    a = _average_ranks(left)
    b = _average_ranks(right)
    if np.std(a) == 0 or np.std(b) == 0:
        raise MetricError("Spearman input is constant")
    return float(np.corrcoef(a, b)[0, 1])


def donor_distance_spearman(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if (
        left.ndim != 2 or right.ndim != 2
        or left.shape[0] != right.shape[0] or len(left) < 3
    ):
        raise MetricError("distance comparison requires matched donor representations")
    upper = np.triu_indices(len(left), k=1)
    return spearman_correlation(
        pairwise_euclidean(left, left)[upper], pairwise_euclidean(right, right)[upper]
    )


def macro_f1(truth: Sequence[str], prediction: Sequence[str]) -> float:
    truth = np.asarray(truth, dtype=object)
    prediction = np.asarray(prediction, dtype=object)
    if truth.shape != prediction.shape or not len(truth):
        raise MetricError("invalid classification vectors")
    labels = sorted(set(truth) | set(prediction))
    scores = []
    for label in labels:
        tp = np.sum((truth == label) & (prediction == label))
        fp = np.sum((truth != label) & (prediction == label))
        fn = np.sum((truth == label) & (prediction != label))
        denominator = 2 * tp + fp + fn
        scores.append(0.0 if denominator == 0 else 2 * tp / denominator)
    return float(np.mean(scores))


def positive_class_f1(truth: Sequence[bool], prediction: Sequence[bool]) -> float:
    """F1 for the prespecified positive lineage, without a dominant negative class."""
    truth = np.asarray(truth, dtype=bool)
    prediction = np.asarray(prediction, dtype=bool)
    if truth.shape != prediction.shape or not len(truth) or not truth.any():
        raise MetricError("positive-class F1 requires matched vectors with positives")
    tp = int(np.sum(truth & prediction))
    fp = int(np.sum(~truth & prediction))
    fn = int(np.sum(truth & ~prediction))
    denominator = 2 * tp + fp + fn
    return 0.0 if denominator == 0 else 2 * tp / denominator


def mean_neighborhood_jaccard(before: np.ndarray, after: np.ndarray, k: int = 30) -> float:
    if before.shape != after.shape or len(before) <= k:
        raise MetricError("neighborhood inputs must be matched with n > k")
    before_d = pairwise_euclidean(before, before)
    after_d = pairwise_euclidean(after, after)
    np.fill_diagonal(before_d, np.inf)
    np.fill_diagonal(after_d, np.inf)
    before_n = np.argpartition(before_d, kth=k - 1, axis=1)[:, :k]
    after_n = np.argpartition(after_d, kth=k - 1, axis=1)[:, :k]
    scores = []
    for left, right in zip(before_n, after_n):
        intersection = len(set(left) & set(right))
        scores.append(intersection / (2 * k - intersection))
    return float(np.mean(scores))


def paired_shift_bootstrap(
    candidate_reference: np.ndarray,
    candidate_control: np.ndarray,
    baseline_reference: np.ndarray,
    baseline_control: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, float]:
    """Paired donor bootstrap of proportional shift-control improvement."""
    if candidate_reference.shape[0] != baseline_reference.shape[0] or candidate_control.shape[0] != baseline_control.shape[0]:
        raise MetricError("candidate and baseline donor rosters must match")
    if min(len(candidate_reference), len(candidate_control)) < 2:
        raise MetricError("bootstrap requires at least two donors per group")
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(replicates):
        r = rng.integers(0, len(candidate_reference), len(candidate_reference))
        c = rng.integers(0, len(candidate_control), len(candidate_control))
        # Duplicated donor draws are valid bootstrap observations; add tiny
        # deterministic jitter only if every within distance collapses to zero.
        try:
            candidate = normalized_shift_control(candidate_reference[r], candidate_control[c])
            baseline = normalized_shift_control(baseline_reference[r], baseline_control[c])
        except MetricError:
            continue
        if baseline > 0:
            values.append((baseline - candidate) / baseline)
    if len(values) < max(100, int(0.9 * replicates)):
        raise MetricError("too many degenerate donor-bootstrap replicates")
    values = np.asarray(values)
    return {
        "estimate": float(np.mean(values)),
        "standard_error": float(np.std(values, ddof=1)),
        "ci_low": float(np.quantile(values, 0.025)),
        "ci_high": float(np.quantile(values, 0.975)),
        "valid_replicates": int(len(values)),
    }


def bh_adjust(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or np.any((values < 0) | (values > 1) | ~np.isfinite(values)):
        raise MetricError("BH inputs must be finite p-values")
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(values)
    result[order] = np.minimum(ranked, 1.0)
    return result
