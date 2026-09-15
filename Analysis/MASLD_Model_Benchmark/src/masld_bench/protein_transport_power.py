"""Prespecified detection floors for the protein_transport task.

Every function here is a pure function of the *design*: sample sizes, class counts,
fold counts, the multiplicity family size and alpha. Nothing in this module can read
a file, a socket, or a subprocess, so no threshold it returns can have been informed
by an observed result. A power figure computed after seeing the effect it is meant to
bound is circular; the structural guarantee is enforced by an AST test rather than
promised in a docstring.

The module is deliberately import-poor. It reuses the repository's own
`average_precision` so the average-precision null is the same estimator the evaluator
will later use, and both that module and `metrics` are themselves stdlib-only and
I/O-free, so the whole import closure stays blind.
"""

from __future__ import annotations

import math
from random import Random
from typing import Mapping, Sequence

from .evaluators.metrics import average_precision

TWO_SIDED = 2


class PowerError(RuntimeError):
    """Raised when a detection floor is requested for a design that cannot support it."""


# Acklam's rational approximation to the inverse standard normal CDF. Used instead of
# scipy so the module keeps a stdlib-only import closure; accurate to about 1e-9,
# far finer than any threshold reported here.
_A = (-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
      1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00)
_B = (-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
      6.680131188771972e01, -1.328068155288572e01)
_C = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
      -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00)
_D = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
      3.754408661907416e00)


def normal_quantile(probability: float) -> float:
    """Inverse standard normal CDF."""
    if not 0.0 < probability < 1.0:
        raise PowerError("a normal quantile needs a probability strictly inside (0, 1)")
    low, high = 0.02425, 1.0 - 0.02425
    if probability < low:
        q = math.sqrt(-2.0 * math.log(probability))
        return (((((_C[0] * q + _C[1]) * q + _C[2]) * q + _C[3]) * q + _C[4]) * q + _C[5]) / (
            (((_D[0] * q + _D[1]) * q + _D[2]) * q + _D[3]) * q + 1.0
        )
    if probability > high:
        q = math.sqrt(-2.0 * math.log(1.0 - probability))
        return -(((((_C[0] * q + _C[1]) * q + _C[2]) * q + _C[3]) * q + _C[4]) * q + _C[5]) / (
            (((_D[0] * q + _D[1]) * q + _D[2]) * q + _D[3]) * q + 1.0
        )
    q = probability - 0.5
    r = q * q
    return (((((_A[0] * r + _A[1]) * r + _A[2]) * r + _A[3]) * r + _A[4]) * r + _A[5]) * q / (
        ((((_B[0] * r + _B[1]) * r + _B[2]) * r + _B[3]) * r + _B[4]) * r + 1.0
    )


def bh_single_discovery_alpha(alpha: float, family_size: int) -> float:
    """The alpha a single test must clear to survive BH as the only discovery.

    Benjamini-Hochberg allows rank k at alpha * k / m. A test that is the sole
    discovery in a family of m is rank 1, so it faces alpha / m. That is the honest
    design bar: it is what a lane must be powered for if it cannot assume its
    neighbours will also fire.
    """
    if family_size < 1:
        raise PowerError("a multiplicity family holds at least one test")
    return alpha / float(family_size)


def detectable_spearman(
    n: int, *, alpha: float = 0.05, power: float = 0.80, tie_free_correction: bool = True
) -> float:
    """Smallest absolute Spearman correlation detectable at the given n and alpha.

    Uses the Fisher z transform with the 1.06 variance inflation appropriate to a rank
    correlation. This is an upper bound on optimism only in the tie-free case; a heavily
    tied outcome compresses the achievable correlation further, so the returned figure
    should be read as the floor under the most favourable tie structure.
    """
    if n < 4:
        raise PowerError("a rank correlation floor needs at least 4 observations")
    inflation = math.sqrt(1.06) if tie_free_correction else 1.0
    standard_error = inflation / math.sqrt(n - 3)
    z = (normal_quantile(1.0 - alpha / TWO_SIDED) + normal_quantile(power)) * standard_error
    return math.tanh(z)


def _auc_null_standard_error(n_positive: int, n_negative: int) -> float:
    return math.sqrt((n_positive + n_negative + 1) / (12.0 * n_positive * n_negative))


def _auc_alternative_standard_error(auc: float, n_positive: int, n_negative: int) -> float:
    """Hanley and McNeil's exponential-approximation variance for AUC."""
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    variance = (
        auc * (1.0 - auc)
        + (n_positive - 1) * (q1 - auc * auc)
        + (n_negative - 1) * (q2 - auc * auc)
    ) / (n_positive * n_negative)
    return math.sqrt(variance)


def detectable_auroc(
    n_positive: int, n_negative: int, *, alpha: float = 0.05, power: float = 0.80
) -> float:
    """Smallest AUROC detectable against 0.5 at the given class counts and alpha."""
    if n_positive < 1 or n_negative < 1:
        raise PowerError("an AUROC floor needs both classes present")
    z_alpha = normal_quantile(1.0 - alpha / TWO_SIDED)
    z_power = normal_quantile(power)
    null_se = _auc_null_standard_error(n_positive, n_negative)
    low, high = 0.5, 0.999999
    for _ in range(200):
        middle = (low + high) / 2.0
        required = z_alpha * null_se + z_power * _auc_alternative_standard_error(
            middle, n_positive, n_negative
        )
        if (middle - 0.5) >= required:
            high = middle
        else:
            low = middle
    return high


def average_precision_null_quantile(
    n_positive: int,
    n_negative: int,
    *,
    quantile: float,
    n_draws: int = 200_000,
    seed: int = 20260825,
) -> float:
    """The operating bar an average precision must clear under a continuous random scorer.

    Class prevalence is NOT this bar. Average precision is upward biased at small n, so
    a scorer that has learned nothing still averages well above prevalence, and scoring
    against prevalence manufactures significance.
    """
    if n_positive < 1 or n_negative < 1:
        raise PowerError("an average precision null needs both classes present")
    if not 0.0 < quantile < 1.0:
        raise PowerError("a quantile lies strictly inside (0, 1)")
    labels = [1] * n_positive + [0] * n_negative
    rng = Random(seed)
    draws = []
    for _ in range(n_draws):
        scores = [rng.random() for _ in labels]
        draws.append(average_precision(labels, scores))
    draws.sort()
    index = min(int(math.ceil(quantile * n_draws)) - 1, n_draws - 1)
    return draws[max(index, 0)]


def average_precision_null_summary(
    n_positive: int,
    n_negative: int,
    *,
    quantiles: Sequence[float],
    n_draws: int = 100_000,
    seed: int = 20260825,
) -> dict:
    """Mean and several quantiles of the average-precision null from ONE simulation.

    Drawing the null once and reading every statistic off the same sorted sample is
    both cheaper than repeated simulation and more coherent: the mean and the bars
    then describe the same realised null rather than three independent ones.
    """
    if n_positive < 1 or n_negative < 1:
        raise PowerError("an average precision null needs both classes present")
    for quantile in quantiles:
        if not 0.0 < quantile < 1.0:
            raise PowerError("a quantile lies strictly inside (0, 1)")
    labels = [1] * n_positive + [0] * n_negative
    rng = Random(seed)
    draws = [average_precision(labels, [rng.random() for _ in labels]) for _ in range(n_draws)]
    total = math.fsum(draws)
    draws.sort()
    quantile_values = {}
    for quantile in quantiles:
        index = min(int(math.ceil(quantile * n_draws)) - 1, n_draws - 1)
        quantile_values[quantile] = draws[max(index, 0)]
    return {
        "mean": total / float(n_draws),
        "quantiles": quantile_values,
        "n_draws": n_draws,
        "seed": seed,
        "n_positive": n_positive,
        "n_negative": n_negative,
        "prevalence": prevalence(n_positive, n_negative),
    }


def average_precision_null_mean(
    n_positive: int, n_negative: int, *, n_draws: int = 200_000, seed: int = 20260825
) -> float:
    """Mean average precision of a continuous random scorer, for contrast with prevalence."""
    labels = [1] * n_positive + [0] * n_negative
    rng = Random(seed)
    total = 0.0
    for _ in range(n_draws):
        total += average_precision(labels, [rng.random() for _ in labels])
    return total / float(n_draws)


def prevalence(n_positive: int, n_negative: int) -> float:
    """Reported only so the gap to the real null can be shown. Never a baseline."""
    return n_positive / float(n_positive + n_negative)


def rank_endpoint_is_estimable(distinct_values_per_fold: Sequence[int]) -> bool:
    """A rank endpoint needs a non-constant held-out vector in every fold, nothing more.

    It does not need every level present in every fold. An absent level costs resolution
    in that part of the range, not definedness.
    """
    return bool(distinct_values_per_fold) and all(
        count >= 2 for count in distinct_values_per_fold
    )


def classification_endpoint_per_fold_estimability(
    level_counts_per_fold: Mapping[str, Sequence[int]],
) -> dict:
    """Which levels vanish from which folds, and whether per-fold scoring survives.

    A level absent from a fold makes a per-fold metric undefined for that level. Pooled
    out-of-fold scoring removes the problem; averaging per-fold values does not.
    """
    if not level_counts_per_fold:
        raise PowerError("estimability needs at least one level")
    empty = {
        level: [index for index, count in enumerate(counts) if count == 0]
        for level, counts in level_counts_per_fold.items()
    }
    vanishing = {level: folds for level, folds in empty.items() if folds}
    return {
        "levels_absent_from_some_fold": sorted(vanishing),
        "empty_level_fold_cells": sum(len(folds) for folds in vanishing.values()),
        "per_fold_metric_estimable": not vanishing,
        "pooled_out_of_fold_estimable": True,
    }


def smallest_estimable_stratum(level_counts: Mapping[str, int], *, outer_folds: int) -> dict:
    """Flag levels too rare to appear in most folds, e.g. Kleiner F3 at 3 of 58.

    A level with fewer members than the fold count cannot appear in every fold even
    under a perfectly balanced split, so any per-fold statement about it is undefined
    by construction rather than by bad luck in the hash.
    """
    if outer_folds < 2:
        raise PowerError("outer folds must be at least 2")
    return {
        level: {
            "count": count,
            "can_appear_in_every_fold": count >= outer_folds,
            "maximum_folds_reachable": min(count, outer_folds),
        }
        for level, count in sorted(level_counts.items())
    }
