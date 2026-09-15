"""
Metric and bootstrap machinery for the encoder benchmark.

NOT written for this release. Every function below is extracted VERBATIM by
AST line range from the repository sources listed here, so the released copy
and the code that produced the published Panel 7F numbers are the same text.
test_release.sh re-derives the published deltas through this copy; if the copy
had drifted, that assertion would fail.

Vendored from:
  scripts/score_cell_baselines_study_50000.py
    sha256 39ce00ce351c096d411a0867a75b59abf77101fb48999e0995d0314cb0b09df0
    build_multiplicities, sufficient_statistics, endpoints_from_multiplicities, interval
  src/masld_bench/evaluators/metrics.py
    sha256 85d048bb8817fc8cb104b3d756981a6cd221415d7095809a07b756608030a078
    donor_class_balanced_weights, weighted_macro_f1

The donor is the unit of resampling. Cells are not independent: 50,000 cells
come from 102 donors, and a cell-level bootstrap would understate every
interval. build_multiplicities is study-stratified and preserves each study's
evaluable-class roster, and the SAME multiplicities are used for every model
so that deltas are paired on the same draws.
"""
from __future__ import annotations

from collections import defaultdict
import math
from typing import Any, Hashable, Mapping, Sequence

import numpy as np


ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
EXPECTED_ROWS = 50_000
EXPECTED_DONORS = 102
EXPECTED_STUDIES = 7
CALIBRATION_BINS = 10
BOOTSTRAP_RESAMPLES = 10_000
BOOTSTRAP_SEED = 20260824


class CellStudyScoreError(ValueError):
    """Raised when frozen predictions cannot support an endpoint."""


class MetricError(ValueError):
    """Raised when standardized predictions cannot support an endpoint."""


def _finite_values(values: Sequence[float], name: str) -> list[float]:
    checked: list[float] = []
    for value in values:
        if isinstance(value, bool):
            raise MetricError(f"{name} must be numeric, not boolean")
        try:
            number = float(value)
        except (TypeError, ValueError) as error:
            raise MetricError(f"{name} must contain numeric values") from error
        if not math.isfinite(number):
            raise MetricError(f"{name} must contain only finite values")
        checked.append(number)
    return checked


def build_multiplicities(
    donor_studies: np.ndarray,
    present: np.ndarray,
    *,
    n_resamples: int,
    seed: int,
) -> np.ndarray:
    """Study-stratified donor bootstrap with fixed evaluable-class coverage."""

    if donor_studies.ndim != 1 or present.shape[0] != len(donor_studies):
        raise CellStudyScoreError("bootstrap donor strata differ")
    rng = np.random.default_rng(seed)
    result = np.zeros((n_resamples, len(donor_studies)), dtype=np.int16)
    for study in sorted(set(donor_studies.tolist())):
        donor_indices = np.flatnonzero(donor_studies == study)
        class_roster = present[donor_indices].any(axis=0)
        remaining = np.arange(n_resamples)
        attempts = 0
        while len(remaining):
            attempts += 1
            if attempts > 10_000:
                raise CellStudyScoreError("class-preserving donor bootstrap stalled")
            draws = rng.integers(
                0, len(donor_indices), size=(len(remaining), len(donor_indices))
            )
            covered = present[donor_indices[draws]].any(axis=1)
            accepted = np.all(covered[:, class_roster], axis=1)
            accepted_rows = remaining[accepted]
            accepted_draws = draws[accepted]
            if len(accepted_rows):
                row_axis = np.repeat(accepted_rows, len(donor_indices))
                column_axis = donor_indices[accepted_draws.reshape(-1)]
                np.add.at(result, (row_axis, column_axis), 1)
            remaining = remaining[~accepted]
    return result

def sufficient_statistics(
    truth: np.ndarray,
    predicted: np.ndarray,
    probabilities: np.ndarray,
    donor_index: np.ndarray,
    *,
    donors: int,
) -> dict[str, np.ndarray]:
    classes = len(ROSTER)
    present = np.zeros((donors, classes), dtype=bool)
    confusion = np.zeros((donors, classes, classes), dtype=np.float64)
    brier = np.zeros((donors, classes), dtype=np.float64)
    bin_mass = np.zeros((donors, classes, CALIBRATION_BINS), dtype=np.float64)
    bin_confidence = np.zeros_like(bin_mass)
    bin_correct = np.zeros_like(bin_mass)
    confidence = probabilities[np.arange(len(truth)), predicted]
    correct = (truth == predicted).astype(np.float64)
    bin_ids = np.minimum(
        (confidence * CALIBRATION_BINS).astype(np.int64), CALIBRATION_BINS - 1
    )
    losses = np.sum(
        (probabilities - np.eye(classes, dtype=np.float64)[truth]) ** 2,
        axis=1,
    )
    for donor in range(donors):
        for truth_class in range(classes):
            selected = (donor_index == donor) & (truth == truth_class)
            count = int(selected.sum())
            if not count:
                continue
            present[donor, truth_class] = True
            confusion[donor, truth_class] = (
                np.bincount(predicted[selected], minlength=classes) / count
            )
            brier[donor, truth_class] = float(losses[selected].mean())
            for bin_id in range(CALIBRATION_BINS):
                in_bin = selected & (bin_ids == bin_id)
                bin_mass[donor, truth_class, bin_id] = in_bin.sum() / count
                bin_confidence[donor, truth_class, bin_id] = confidence[in_bin].sum() / count
                bin_correct[donor, truth_class, bin_id] = correct[in_bin].sum() / count
    return {
        "present": present,
        "confusion": confusion,
        "brier": brier,
        "bin_mass": bin_mass,
        "bin_confidence": bin_confidence,
        "bin_correct": bin_correct,
    }

def endpoints_from_multiplicities(
    stats: Mapping[str, np.ndarray],
    multiplicities: np.ndarray,
    donor_indices: np.ndarray,
    class_roster: np.ndarray,
) -> dict[str, np.ndarray]:
    mult = multiplicities[:, donor_indices].astype(np.float64)
    present = stats["present"][donor_indices][:, class_roster]
    denominators = mult @ present.astype(np.float64)
    if np.any(denominators <= 0.0):
        raise CellStudyScoreError("bootstrap lost a frozen evaluation class")
    confusion_sum = np.einsum(
        "bd,dcp->bcp",
        mult,
        stats["confusion"][donor_indices][:, class_roster, :],
        optimize=True,
    )
    weighted_confusion = (
        confusion_sum / denominators[:, :, np.newaxis] / int(class_roster.sum())
    )
    target_positions = np.flatnonzero(class_roster)
    true_positive = np.stack(
        [weighted_confusion[:, position, label] for position, label in enumerate(target_positions)],
        axis=1,
    )
    false_positive = np.stack(
        [weighted_confusion[:, :, label].sum(axis=1) for label in target_positions],
        axis=1,
    ) - true_positive
    false_negative = weighted_confusion.sum(axis=2) - true_positive
    denominator = 2.0 * true_positive + false_positive + false_negative
    class_f1 = np.divide(
        2.0 * true_positive,
        denominator,
        out=np.zeros_like(true_positive),
        where=denominator > 0.0,
    )
    brier_sum = mult @ stats["brier"][donor_indices][:, class_roster]
    brier = np.mean(brier_sum / denominators, axis=1)
    confidence_sum = np.einsum(
        "bd,dck->bck",
        mult,
        stats["bin_confidence"][donor_indices][:, class_roster, :],
        optimize=True,
    )
    correct_sum = np.einsum(
        "bd,dck->bck",
        mult,
        stats["bin_correct"][donor_indices][:, class_roster, :],
        optimize=True,
    )
    weighted_confidence = (
        confidence_sum / denominators[:, :, np.newaxis] / int(class_roster.sum())
    ).sum(axis=1)
    weighted_correct = (
        correct_sum / denominators[:, :, np.newaxis] / int(class_roster.sum())
    ).sum(axis=1)
    ece = np.abs(weighted_correct - weighted_confidence).sum(axis=1)
    return {
        "macro_f1": class_f1.mean(axis=1),
        "class_f1": class_f1,
        "brier": brier,
        "top_label_ece": ece,
    }

def interval(values: np.ndarray) -> dict[str, Any]:
    if values.shape != (BOOTSTRAP_RESAMPLES,) or not np.all(np.isfinite(values)):
        raise CellStudyScoreError("bootstrap distribution differs")
    return {
        "lower": float(np.quantile(values, 0.025, method="linear")),
        "median": float(np.quantile(values, 0.5, method="linear")),
        "upper": float(np.quantile(values, 0.975, method="linear")),
        "n_resamples": BOOTSTRAP_RESAMPLES,
        "confidence_level": 0.95,
        "seed": BOOTSTRAP_SEED,
    }

def donor_class_balanced_weights(
    donor_ids: Sequence[Hashable], labels: Sequence[Hashable]
) -> list[float]:
    """Give every class equal mass, then every donor within that class equal mass.

    This is deliberately not a cell-balanced weight and not a flat
    donor-by-class-stratum weight.  The latter overweights classes observed in
    more donors.  Here, each class contributes ``1 / n_classes``; within a
    class, every donor containing that class contributes equally; and cells
    only divide their donor-by-class share.
    """

    if len(donor_ids) != len(labels) or not labels:
        raise MetricError("donor_ids and labels must be non-empty and aligned")
    counts: dict[tuple[Hashable, Hashable], int] = defaultdict(int)
    donors_by_class: dict[Hashable, set[Hashable]] = defaultdict(set)
    for donor_id, label in zip(donor_ids, labels, strict=True):
        counts[(donor_id, label)] += 1
        donors_by_class[label].add(donor_id)
    class_count = len(donors_by_class)
    return [
        1.0
        / (
            class_count
            * len(donors_by_class[label])
            * counts[(donor_id, label)]
        )
        for donor_id, label in zip(donor_ids, labels, strict=True)
    ]

def weighted_macro_f1(
    observed: Sequence[Hashable],
    predicted: Sequence[Hashable],
    weights: Sequence[float],
) -> tuple[float, dict[Hashable, float]]:
    if not (len(observed) == len(predicted) == len(weights)) or not observed:
        raise MetricError("classification arrays must be non-empty and aligned")
    checked_weights = _finite_values(weights, "classification weights")
    if any(weight < 0.0 for weight in checked_weights) or sum(checked_weights) <= 0.0:
        raise MetricError("classification weights must be nonnegative with positive mass")
    classes = sorted(set(observed), key=str)
    class_scores: dict[Hashable, float] = {}
    for label in classes:
        true_positive = sum(
            weight
            for truth, guess, weight in zip(
                observed, predicted, checked_weights, strict=True
            )
            if truth == label and guess == label
        )
        false_positive = sum(
            weight
            for truth, guess, weight in zip(
                observed, predicted, checked_weights, strict=True
            )
            if truth != label and guess == label
        )
        false_negative = sum(
            weight
            for truth, guess, weight in zip(
                observed, predicted, checked_weights, strict=True
            )
            if truth == label and guess != label
        )
        denominator = 2 * true_positive + false_positive + false_negative
        class_scores[label] = 0.0 if denominator == 0 else 2 * true_positive / denominator
    return sum(class_scores.values()) / len(class_scores), class_scores
