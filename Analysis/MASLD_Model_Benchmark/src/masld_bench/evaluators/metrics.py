"""Dependency-free benchmark metrics computed from frozen predictions."""

from __future__ import annotations

from collections import defaultdict
import math
from typing import Any, Hashable, Iterable, Mapping, Sequence


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


def multiclass_brier_score(
    observed: Sequence[Hashable],
    probabilities: Sequence[Mapping[Hashable, float]],
    weights: Sequence[float],
) -> float:
    if not (len(observed) == len(probabilities) == len(weights)) or not observed:
        raise MetricError("Brier inputs must be non-empty and aligned")
    checked_weights = _finite_values(weights, "classification weights")
    if any(weight < 0.0 for weight in checked_weights):
        raise MetricError("classification weights must be nonnegative")
    classes = set(observed)
    checked_probabilities: list[dict[Hashable, float]] = []
    for row in probabilities:
        if not isinstance(row, Mapping) or not row:
            raise MetricError("each probability row must be a non-empty mapping")
        classes.update(row)
        values = _finite_values(tuple(row.values()), "class probabilities")
        total = sum(values)
        if abs(total - 1.0) > 1e-6 or any(value < 0.0 for value in values):
            raise MetricError("each probability row must be nonnegative and sum to one")
        checked_probabilities.append(dict(zip(row, values, strict=True)))
    total_weight = sum(checked_weights)
    if total_weight <= 0:
        raise MetricError("classification weights must sum to a positive value")
    loss = 0.0
    for truth, row, weight in zip(
        observed, checked_probabilities, checked_weights, strict=True
    ):
        loss += weight * sum(
            (float(row.get(label, 0.0)) - float(label == truth)) ** 2 for label in classes
        )
    return loss / total_weight


def _average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        average = ((start + 1) + stop) / 2.0
        for offset in range(start, stop):
            ranks[order[offset]] = average
        start = stop
    return ranks


def pearson_correlation(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        raise MetricError("correlation requires aligned vectors with at least two values")
    left_values = _finite_values(left, "left correlation vector")
    right_values = _finite_values(right, "right correlation vector")
    left_mean = sum(left_values) / len(left_values)
    right_mean = sum(right_values) / len(right_values)
    numerator = sum(
        (x - left_mean) * (y - right_mean)
        for x, y in zip(left_values, right_values, strict=True)
    )
    left_ss = sum((x - left_mean) ** 2 for x in left_values)
    right_ss = sum((y - right_mean) ** 2 for y in right_values)
    if left_ss == 0 or right_ss == 0:
        raise MetricError("correlation is undefined for a constant vector")
    return numerator / math.sqrt(left_ss * right_ss)


def spearman_correlation(left: Sequence[float], right: Sequence[float]) -> float:
    return pearson_correlation(_average_ranks(left), _average_ranks(right))


def fisher_z_mean(correlations: Sequence[float]) -> float:
    if not correlations:
        raise MetricError("at least one correlation is required")
    checked = _finite_values(correlations, "correlations")
    if any(value < -1.0 or value > 1.0 for value in checked):
        raise MetricError("correlations must lie within [-1, 1]")
    epsilon = 1e-12
    z_values = [
        math.atanh(max(-1 + epsilon, min(1 - epsilon, value)))
        for value in checked
    ]
    return math.tanh(sum(z_values) / len(z_values))


def average_precision(observed_binary: Sequence[int], scores: Sequence[float]) -> float:
    """Tie-stable non-interpolated average precision."""

    if len(observed_binary) != len(scores) or not scores:
        raise MetricError("AUPRC inputs must be non-empty and aligned")
    if any(isinstance(value, bool) or value not in {0, 1} for value in observed_binary):
        raise MetricError("AUPRC labels must be binary")
    checked_scores = _finite_values(scores, "AUPRC scores")
    positives = sum(observed_binary)
    if positives == 0:
        raise MetricError("AUPRC is undefined without a positive observation")
    grouped: dict[float, list[int]] = defaultdict(list)
    for truth, score in zip(observed_binary, checked_scores, strict=True):
        grouped[score].append(truth)
    true_positive = 0
    false_positive = 0
    previous_recall = 0.0
    area = 0.0
    for score in sorted(grouped, reverse=True):
        group = grouped[score]
        true_positive += sum(group)
        false_positive += len(group) - sum(group)
        recall = true_positive / positives
        precision = true_positive / (true_positive + false_positive)
        area += (recall - previous_recall) * precision
        previous_recall = recall
    return area


def paired_cosine_retrieval_accuracy(
    query_vectors: Sequence[Sequence[float]],
    target_vectors: Sequence[Sequence[float]],
    query_pair_ids: Sequence[Hashable],
    target_pair_ids: Sequence[Hashable],
    *,
    top_k: int = 1,
) -> float:
    """Measure same-entity cross-assay retrieval without treating pairs as replicates."""

    if not query_vectors or not target_vectors:
        raise MetricError("retrieval vectors must be non-empty")
    if len(query_vectors) != len(query_pair_ids) or len(target_vectors) != len(
        target_pair_ids
    ):
        raise MetricError("retrieval vectors and pair identifiers must align")
    if len(set(query_pair_ids)) != len(query_pair_ids) or len(set(target_pair_ids)) != len(
        target_pair_ids
    ):
        raise MetricError("same-nucleus retrieval pair identifiers must be unique")
    if set(query_pair_ids) != set(target_pair_ids):
        raise MetricError("query and target retrieval pairs must match exactly")
    if isinstance(top_k, bool) or not isinstance(top_k, int) or not 1 <= top_k <= len(
        target_vectors
    ):
        raise MetricError("top_k must be between one and the target count")
    dimensions = {len(vector) for vector in (*query_vectors, *target_vectors)}
    if len(dimensions) != 1 or next(iter(dimensions)) < 1:
        raise MetricError("retrieval embeddings must share a positive dimension")

    def normalized(vector: Sequence[float]) -> tuple[float, ...]:
        values = tuple(float(value) for value in vector)
        if any(not math.isfinite(value) for value in values):
            raise MetricError("retrieval embeddings must be finite")
        norm = math.sqrt(sum(value * value for value in values))
        if norm <= 0:
            raise MetricError("retrieval embeddings may not contain a zero vector")
        return tuple(value / norm for value in values)

    normalized_queries = tuple(normalized(vector) for vector in query_vectors)
    normalized_targets = tuple(normalized(vector) for vector in target_vectors)
    hits = 0
    for query, pair_id in zip(normalized_queries, query_pair_ids, strict=True):
        ranked = sorted(
            zip(target_pair_ids, normalized_targets, strict=True),
            key=lambda item: (
                -sum(left * right for left, right in zip(query, item[1], strict=True)),
                str(item[0]),
            ),
        )
        hits += pair_id in {candidate_id for candidate_id, _ in ranked[:top_k]}
    return hits / len(query_pair_ids)


def multinomial_deviance(observed: Sequence[float], predicted: Sequence[float]) -> float:
    if len(observed) != len(predicted) or not observed:
        raise MetricError("profile vectors must be non-empty and aligned")
    observed_values = _finite_values(observed, "observed profile")
    predicted_values = _finite_values(predicted, "predicted profile")
    if any(value < 0 for value in observed_values) or any(
        value < 0 for value in predicted_values
    ):
        raise MetricError("profile values must be nonnegative")
    observed_total = sum(observed_values)
    predicted_total = sum(predicted_values)
    if observed_total <= 0 or predicted_total <= 0:
        raise MetricError("profile vectors must have positive totals")
    expected = [observed_total * value / predicted_total for value in predicted_values]
    if any(
        value <= 0 and truth > 0
        for truth, value in zip(observed_values, expected, strict=True)
    ):
        return math.inf
    return 2.0 * sum(
        truth * math.log(truth / value)
        for truth, value in zip(observed_values, expected, strict=True)
        if truth > 0
    )


def relative_deviance_reduction(
    observed_profiles: Sequence[Sequence[float]],
    model_profiles: Sequence[Sequence[float]],
    baseline_profiles: Sequence[Sequence[float]],
) -> float:
    if not (
        len(observed_profiles) == len(model_profiles) == len(baseline_profiles)
    ) or not observed_profiles:
        raise MetricError("profile collections must be non-empty and aligned")
    model_deviance = sum(
        multinomial_deviance(observed, predicted)
        for observed, predicted in zip(observed_profiles, model_profiles, strict=True)
    )
    baseline_deviance = sum(
        multinomial_deviance(observed, predicted)
        for observed, predicted in zip(observed_profiles, baseline_profiles, strict=True)
    )
    if not math.isfinite(model_deviance):
        raise MetricError("model deviance must be finite")
    if baseline_deviance <= 0 or not math.isfinite(baseline_deviance):
        raise MetricError("baseline deviance must be finite and positive")
    return (baseline_deviance - model_deviance) / baseline_deviance
