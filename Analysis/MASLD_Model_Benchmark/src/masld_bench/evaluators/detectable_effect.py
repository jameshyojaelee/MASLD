"""Detectable-effect curves for pooled ranking metrics at a fixed design.

A power statement earns its place only if it is derived from the design rather
than from the result. Everything here is a function of the frozen design alone:
the outcome vector and its denominators, the decision rule, and the multiplicity
family. **No observed model performance enters any calculation in this module.**
An effect is planted at a stated size and the detection rate is measured; the
realised effect of any candidate is never read.

Two curves answer two different questions.

``information_ceiling``
    Give an *oracle* predictor whose true association with the outcome is known,
    and ask how often the sample metric clears the decision threshold. This is
    the ceiling: no model, however good, can beat it at this n, because the limit
    is the sampling noise of the metric and not the quality of the learner. A
    design that cannot detect an effect here cannot detect it at all.

``realised_power``
    The same question after cross-fitted estimation noise is added, which is what
    an actual model faces. It is always below the ceiling, and the gap is the
    price of having to learn the predictor rather than being handed it.

Effect sizes are reported on the scale a reader can interpret — the mean sample
Spearman or the mean sample average precision the planted effect actually
produces — rather than on the latent parameter used to generate it.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, sqrt
from random import Random
from typing import Sequence

from .metrics import MetricError, average_precision, spearman_correlation


#: Decision rules the campaign actually operates under.
#:
#: ``nominal`` ignores multiplicity and is an upper bound on power.
#: ``bh_single_true`` is Benjamini-Hochberg at q0.05 within the twelve-test
#: confirmatory family when exactly one test is truly non-null. BH rejects the
#: smallest p only if it clears alpha/m, so this reduces to Bonferroni and is a
#: lower bound on power. If more scopes were truly non-null BH would be less
#: stringent, so the two rules bracket the operating characteristic rather than
#: pretending to a single number.
NOMINAL_ALPHA = 0.05
CONFIRMATORY_FAMILY_SIZE = 12
BH_SINGLE_TRUE_ALPHA = NOMINAL_ALPHA / CONFIRMATORY_FAMILY_SIZE

TARGET_POWER = 0.80
DEFAULT_NULL_DRAWS = 40_000
DEFAULT_REPLICATES = 4_000
DEFAULT_SEED = 20260825


@dataclass(frozen=True)
class PowerPoint:
    """Power at one planted effect size under both decision rules."""

    latent_parameter: float
    mean_sample_effect: float
    secondary_effect: float
    power_nominal: float
    power_bh_single_true: float
    replicates: int
    n: int

    def to_dict(self) -> dict[str, object]:
        return {
            "latent_parameter": self.latent_parameter,
            "mean_sample_effect": self.mean_sample_effect,
            "secondary_effect": self.secondary_effect,
            "power_nominal": self.power_nominal,
            "power_bh_single_true": self.power_bh_single_true,
            "replicates": self.replicates,
            "n": self.n,
        }


def _quantile(ordered: Sequence[float], probability: float) -> float:
    if not ordered:
        raise MetricError("no draws to take a quantile from")
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
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


def _standardised_ranks(outcome: Sequence[float]) -> list[float]:
    ranks = _average_ranks(outcome)
    mean = sum(ranks) / len(ranks)
    spread = sqrt(sum((value - mean) ** 2 for value in ranks) / len(ranks))
    if spread == 0.0:
        raise MetricError("outcome vector is constant")
    return [(value - mean) / spread for value in ranks]


def spearman_critical_values(
    outcome: Sequence[float],
    *,
    n_draws: int = DEFAULT_NULL_DRAWS,
    seed: int = DEFAULT_SEED,
) -> dict[str, float]:
    """Critical |rho| a continuous random scorer clears at each rule's level.

    Computed against the real outcome vector, so its tie structure is respected
    even though ties in the outcome barely move a rank null.
    """

    values = list(outcome)
    rng = Random(seed)
    draws = [
        abs(spearman_correlation([rng.random() for _ in values], values))
        for _ in range(n_draws)
    ]
    draws.sort()
    return {
        "nominal": _quantile(draws, 1.0 - NOMINAL_ALPHA),
        "bh_single_true": _quantile(draws, 1.0 - BH_SINGLE_TRUE_ALPHA),
    }


def auprc_critical_values(
    n_positive: int,
    n_negative: int,
    *,
    n_draws: int = DEFAULT_NULL_DRAWS,
    seed: int = DEFAULT_SEED,
    labels: Sequence[int] | None = None,
) -> dict[str, float]:
    """Average precision a continuous random scorer clears at each rule's level.

    ``labels`` accepts the label vector in its real order. With continuous
    scores the population quantity depends only on the label multiset, but a
    *fixed seed* pairs its draws with label positions, so ordering is a free
    parameter of the estimate at any finite draw count. It is invisible in a
    provenance record that lists only ``n_draws`` and ``seed``, and at ten
    thousand draws it moves the fibrosis p95 by about 0.003. Pass the real
    vector and record which ordering was used.
    """

    if labels is not None:
        labels = list(labels)
        if sum(labels) != n_positive or len(labels) - sum(labels) != n_negative:
            raise MetricError("supplied labels disagree with the stated counts")
    else:
        labels = [1] * n_positive + [0] * n_negative
    rng = Random(seed)
    draws = [
        average_precision(labels, [rng.random() for _ in labels])
        for _ in range(n_draws)
    ]
    draws.sort()
    return {
        "nominal": _quantile(draws, 1.0 - NOMINAL_ALPHA),
        "bh_single_true": _quantile(draws, 1.0 - BH_SINGLE_TRUE_ALPHA),
    }


def spearman_power_curve(
    outcome: Sequence[float],
    latent_grid: Sequence[float],
    *,
    replicates: int = DEFAULT_REPLICATES,
    null_draws: int = DEFAULT_NULL_DRAWS,
    seed: int = DEFAULT_SEED,
    attenuation: float = 1.0,
) -> list[PowerPoint]:
    """Power to detect a planted rank association, as a function of its size.

    ``attenuation`` multiplies the latent signal before noise is added, which is
    how cross-fitted estimation loss is represented without simulating a
    learner. Leave it at one for the information ceiling.
    """

    values = list(outcome)
    critical = spearman_critical_values(values, n_draws=null_draws, seed=seed)
    reference = _standardised_ranks(values)
    rng = Random(seed + 1)
    points: list[PowerPoint] = []
    for latent in latent_grid:
        effective = latent * attenuation
        residual = sqrt(max(0.0, 1.0 - effective * effective))
        observed: list[float] = []
        hits_nominal = 0
        hits_bh = 0
        for _ in range(replicates):
            predictor = [
                effective * anchor + residual * rng.gauss(0.0, 1.0)
                for anchor in reference
            ]
            sample = spearman_correlation(predictor, values)
            observed.append(sample)
            if abs(sample) >= critical["nominal"]:
                hits_nominal += 1
            if abs(sample) >= critical["bh_single_true"]:
                hits_bh += 1
        points.append(
            PowerPoint(
                latent_parameter=latent,
                mean_sample_effect=sum(observed) / len(observed),
                secondary_effect=float("nan"),
                power_nominal=hits_nominal / replicates,
                power_bh_single_true=hits_bh / replicates,
                replicates=replicates,
                n=len(values),
            )
        )
    return points


def _auroc_from_separation(separation: float) -> float:
    """Two-sample AUROC implied by a unit-variance mean shift."""

    return 0.5 * (1.0 + erf(separation / 2.0))


def auprc_power_curve(
    n_positive: int,
    n_negative: int,
    separation_grid: Sequence[float],
    *,
    replicates: int = DEFAULT_REPLICATES,
    null_draws: int = DEFAULT_NULL_DRAWS,
    seed: int = DEFAULT_SEED,
    attenuation: float = 1.0,
) -> list[PowerPoint]:
    """Power to detect a planted class separation, as a function of its size.

    ``separation_grid`` is in standard-deviation units, so the implied true
    AUROC is ``Phi(d / sqrt(2))`` and is reported alongside.
    """

    labels = [1] * n_positive + [0] * n_negative
    critical = auprc_critical_values(
        n_positive, n_negative, n_draws=null_draws, seed=seed
    )
    rng = Random(seed + 1)
    points: list[PowerPoint] = []
    for separation in separation_grid:
        effective = separation * attenuation
        observed: list[float] = []
        hits_nominal = 0
        hits_bh = 0
        for _ in range(replicates):
            scores = [
                rng.gauss(effective if label else 0.0, 1.0) for label in labels
            ]
            sample = average_precision(labels, scores)
            observed.append(sample)
            if sample >= critical["nominal"]:
                hits_nominal += 1
            if sample >= critical["bh_single_true"]:
                hits_bh += 1
        points.append(
            PowerPoint(
                latent_parameter=separation,
                mean_sample_effect=sum(observed) / len(observed),
                secondary_effect=_auroc_from_separation(effective),
                power_nominal=hits_nominal / replicates,
                power_bh_single_true=hits_bh / replicates,
                replicates=replicates,
                n=len(labels),
            )
        )
    return points


def minimum_detectable_effect(
    points: Sequence[PowerPoint], *, rule: str, target: float = TARGET_POWER
) -> dict[str, float] | None:
    """Interpolate the smallest planted effect reaching the target power."""

    field = "power_nominal" if rule == "nominal" else "power_bh_single_true"
    ordered = sorted(points, key=lambda point: point.latent_parameter)
    previous: PowerPoint | None = None
    for point in ordered:
        power = getattr(point, field)
        if power >= target:
            if previous is None:
                return {
                    "latent_parameter": point.latent_parameter,
                    "sample_effect": point.mean_sample_effect,
                    "secondary_effect": point.secondary_effect,
                    "interpolated": 0.0,
                }
            span = getattr(point, field) - getattr(previous, field)
            weight = 0.0 if span == 0 else (target - getattr(previous, field)) / span
            return {
                "latent_parameter": previous.latent_parameter
                + weight * (point.latent_parameter - previous.latent_parameter),
                "sample_effect": previous.mean_sample_effect
                + weight * (point.mean_sample_effect - previous.mean_sample_effect),
                "secondary_effect": previous.secondary_effect
                + weight * (point.secondary_effect - previous.secondary_effect),
                "interpolated": 1.0,
            }
        previous = point
    return None


def resample_outcome(
    outcome: Sequence[float], size: int, *, seed: int = DEFAULT_SEED
) -> list[float]:
    """Extend an outcome vector to a larger n, preserving its shape and ties.

    Used only for the required-sample-size sweep, so that a hypothetical larger
    cohort keeps the distributional awkwardness of the real one rather than
    being quietly replaced by a well-behaved one.
    """

    values = list(outcome)
    rng = Random(seed)
    return [values[rng.randrange(len(values))] for _ in range(size)]


def required_sample_size(
    outcome: Sequence[float],
    sample_sizes: Sequence[int],
    latent: float,
    *,
    rule: str = "bh_single_true",
    replicates: int = 1_500,
    null_draws: int = 12_000,
    seed: int = DEFAULT_SEED,
    target: float = TARGET_POWER,
) -> list[dict[str, object]]:
    """Power at a fixed planted effect across candidate cohort sizes."""

    field = "power_nominal" if rule == "nominal" else "power_bh_single_true"
    rows: list[dict[str, object]] = []
    for size in sample_sizes:
        extended = resample_outcome(outcome, size, seed=seed + size)
        point = spearman_power_curve(
            extended,
            [latent],
            replicates=replicates,
            null_draws=null_draws,
            seed=seed + size,
        )[0]
        rows.append(
            {
                "n": size,
                "latent_parameter": latent,
                "mean_sample_effect": point.mean_sample_effect,
                "power_nominal": point.power_nominal,
                "power_bh_single_true": point.power_bh_single_true,
                "reaches_target": getattr(point, field) >= target,
            }
        )
    return rows


def required_sample_size_binary(
    prevalence: float,
    sample_sizes: Sequence[int],
    separation: float,
    *,
    rule: str = "bh_single_true",
    replicates: int = 1_500,
    null_draws: int = 12_000,
    seed: int = DEFAULT_SEED,
    target: float = TARGET_POWER,
) -> list[dict[str, object]]:
    """Power at a fixed planted separation across candidate cohort sizes."""

    field = "power_nominal" if rule == "nominal" else "power_bh_single_true"
    rows: list[dict[str, object]] = []
    for size in sample_sizes:
        positives = max(1, round(size * prevalence))
        negatives = size - positives
        point = auprc_power_curve(
            positives,
            negatives,
            [separation],
            replicates=replicates,
            null_draws=null_draws,
            seed=seed + size,
        )[0]
        rows.append(
            {
                "n": size,
                "n_positive": positives,
                "n_negative": negatives,
                "separation": separation,
                "implied_true_auroc": _auroc_from_separation(separation),
                "mean_sample_auprc": point.mean_sample_effect,
                "power_nominal": point.power_nominal,
                "power_bh_single_true": point.power_bh_single_true,
                "reaches_target": getattr(point, field) >= target,
            }
        )
    return rows


#: Boundaries for the plain-language verdict. A power statement that can only
#: ever conclude "underpowered" is not a power statement, so the adequate branch
#: is reachable and is exercised by a unit test.
ADEQUATE_UPPER_BOUND = 60
MODERATE_UPPER_BOUND = 200


def power_verdict(
    smallest_sufficient_n: Sequence[int | None],
    *,
    cohort_n: int,
) -> str:
    """Classify a design from the sample sizes its reference effects need.

    ``smallest_sufficient_n`` holds, for each reference effect and rule, the
    smallest cohort reaching the target power, or ``None`` where no tested size
    reached it.
    """

    reached = [value for value in smallest_sufficient_n if value is not None]
    if not reached:
        return "no_reference_effect_reaches_target_power_within_the_tested_range"
    largest = max(reached)
    if largest <= max(ADEQUATE_UPPER_BOUND, cohort_n):
        return "adequately_powered_at_the_realised_cohort_size"
    if largest < MODERATE_UPPER_BOUND:
        return "underpowered_by_roughly_a_factor_of_two_to_five"
    return "underpowered_by_an_order_of_magnitude"
