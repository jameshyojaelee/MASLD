"""Independent, deterministic perturbation baselines on biological-unit profiles.

This module implements the project formulas documented in the perturbation-native
baseline audit. It contains no copied Systema source code.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence


Vector = tuple[float, ...]
Matrix = tuple[Vector, ...]


class PerturbationBaselineError(ValueError):
    """Raised when a baseline request violates topology or missingness rules."""


@dataclass(frozen=True, slots=True)
class CellCountRow:
    biological_unit: str
    target: str
    guide_family: str
    context: str
    is_control: bool
    counts: Vector


@dataclass(frozen=True, slots=True)
class PseudobulkRow:
    biological_unit: str
    target: str
    context: str
    is_control: bool
    counts: Vector
    cells: int
    guide_families: tuple[str, ...]


def _check_vectors(vectors: Sequence[Sequence[float]], label: str) -> int:
    if not vectors:
        raise PerturbationBaselineError(f"{label} is empty")
    width = len(vectors[0])
    if width == 0 or any(len(row) != width for row in vectors):
        raise PerturbationBaselineError(f"{label} has inconsistent width")
    return width


def mean_vector(vectors: Sequence[Sequence[float]]) -> Vector:
    width = _check_vectors(vectors, "vectors")
    return tuple(sum(float(row[index]) for row in vectors) / len(vectors) for index in range(width))


def pseudobulk_counts(rows: Iterable[CellCountRow]) -> tuple[PseudobulkRow, ...]:
    """Sum cells within biological unit, target, context, and control state."""

    groups: dict[tuple[str, str, str, bool], list[CellCountRow]] = defaultdict(list)
    for row in rows:
        if not row.biological_unit or not row.target or not row.guide_family or not row.context:
            raise PerturbationBaselineError("cell row has an empty topology identifier")
        groups[(row.biological_unit, row.target, row.context, row.is_control)].append(row)
    if not groups:
        raise PerturbationBaselineError("no cell rows")
    output: list[PseudobulkRow] = []
    for key in sorted(groups):
        values = groups[key]
        width = _check_vectors([row.counts for row in values], "cell counts")
        summed = tuple(sum(float(row.counts[index]) for row in values) for index in range(width))
        output.append(
            PseudobulkRow(
                biological_unit=key[0],
                target=key[1],
                context=key[2],
                is_control=key[3],
                counts=summed,
                cells=len(values),
                guide_families=tuple(sorted({row.guide_family for row in values})),
            )
        )
    return tuple(output)


def control_mean(rows: Sequence[PseudobulkRow], context: str) -> Vector:
    """Return the equal-biological-unit mean of control pseudobulks."""

    by_unit: dict[str, list[Vector]] = defaultdict(list)
    for row in rows:
        if row.is_control and row.context == context:
            by_unit[row.biological_unit].append(row.counts)
    if not by_unit:
        raise PerturbationBaselineError("matched training controls are structurally missing")
    unit_means = [mean_vector(by_unit[unit]) for unit in sorted(by_unit)]
    return mean_vector(unit_means)


def perturbed_mean(rows: Sequence[PseudobulkRow], context: str) -> Vector:
    """Return a population perturbation prior with total weight one per unit."""

    by_unit: dict[str, list[Vector]] = defaultdict(list)
    for row in rows:
        if not row.is_control and row.context == context:
            by_unit[row.biological_unit].append(row.counts)
    if not by_unit:
        raise PerturbationBaselineError("training perturbations are structurally missing")
    unit_means = [mean_vector(by_unit[unit]) for unit in sorted(by_unit)]
    return mean_vector(unit_means)


def additive_features(context: Sequence[float], target: Sequence[float]) -> Vector:
    return (1.0, *(float(value) for value in context), *(float(value) for value in target))


def bilinear_features(context: Sequence[float], target: Sequence[float]) -> Vector:
    interaction = tuple(float(left) * float(right) for left in context for right in target)
    return (*additive_features(context, target), *interaction)


def _solve(matrix: list[list[float]], right: list[list[float]]) -> Matrix:
    size = len(matrix)
    width = len(right[0])
    augmented = [matrix[index][:] + right[index][:] for index in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-12:
            raise PerturbationBaselineError("ridge normal equation is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    augmented[row][index] - factor * augmented[column][index]
                    for index in range(size + width)
                ]
    return tuple(tuple(row[size:]) for row in augmented)


def fit_ridge(
    design: Sequence[Sequence[float]],
    responses: Sequence[Sequence[float]],
    *,
    penalty: float,
    weights: Sequence[float] | None = None,
) -> Matrix:
    """Fit deterministic multioutput ridge, leaving the intercept unpenalized."""

    features = _check_vectors(design, "design")
    outputs = _check_vectors(responses, "responses")
    if len(design) != len(responses) or penalty < 0:
        raise PerturbationBaselineError("ridge row count or penalty differs")
    observed_weights = tuple(1.0 for _ in design) if weights is None else tuple(float(value) for value in weights)
    if len(observed_weights) != len(design) or any(value <= 0 for value in observed_weights):
        raise PerturbationBaselineError("ridge weights must be positive and row aligned")
    gram = [[0.0 for _ in range(features)] for _ in range(features)]
    cross = [[0.0 for _ in range(outputs)] for _ in range(features)]
    for row, response, weight in zip(design, responses, observed_weights, strict=True):
        for left in range(features):
            for right_index in range(features):
                gram[left][right_index] += weight * float(row[left]) * float(row[right_index])
            for output in range(outputs):
                cross[left][output] += weight * float(row[left]) * float(response[output])
    for index in range(1, features):
        gram[index][index] += penalty
    return _solve(gram, cross)


def predict_linear(design: Sequence[Sequence[float]], coefficients: Matrix) -> Matrix:
    features = len(coefficients)
    outputs = len(coefficients[0]) if coefficients else 0
    if not coefficients or any(len(row) != outputs for row in coefficients):
        raise PerturbationBaselineError("coefficient matrix differs")
    predictions: list[Vector] = []
    for row in design:
        if len(row) != features:
            raise PerturbationBaselineError("prediction design width differs")
        predictions.append(
            tuple(sum(float(row[index]) * coefficients[index][output] for index in range(features)) for output in range(outputs))
        )
    return tuple(predictions)


def systema_matching_mean(
    components: Sequence[str],
    singleton_training_means: Mapping[str, Sequence[float]],
    population_perturbed_mean: Sequence[float],
    *,
    held_single_component: bool = False,
) -> Vector:
    """Implement the published matching mean after biological-unit aggregation."""

    fallback = tuple(float(value) for value in population_perturbed_mean)
    if not components:
        raise PerturbationBaselineError("Systema components are empty")
    if held_single_component and len(components) == 1:
        return fallback
    matched = [tuple(float(value) for value in singleton_training_means.get(component, fallback)) for component in components]
    return mean_vector(matched)


def strict_training_reference(profile: Sequence[float], training_centroid: Sequence[float]) -> Vector:
    if len(profile) != len(training_centroid) or not profile:
        raise PerturbationBaselineError("strict reference width differs")
    return tuple(float(value) - float(center) for value, center in zip(profile, training_centroid, strict=True))
