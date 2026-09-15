"""Null references for ranking metrics at small sample sizes.

Two facts about average precision decide whether a reported AUPRC means
anything, and neither is visible from the number itself.

**Average precision is upward biased.**  The prevalence of the positive class
is often quoted as the "baseline" a classifier must beat.  It is not the
expected value of average precision under a random scorer.  At small sample
sizes a random scorer averages materially above prevalence, so scoring against
the prevalence line converts noise into apparent skill.

**The permutation null inherits the tie structure of the scores it permutes.**
A model emitting few distinct values cannot spread its permuted average
precision as widely as a continuous scorer can, so its null is tighter.  A
single fixed null is therefore wrong for both: it is too tight for the
continuous case and too loose for the tied one.

The consequence for practice is that a reported AUPRC needs two references,
not one:

``random_score_reference``
    what a *continuous* random scorer reaches on this label vector.  This is
    the honest detection floor, and it belongs beside every reported number.
``permutation_reference``
    what *this model's own score vector* reaches when the labels are shuffled.
    This respects the model's ties and gives the p-value.

This module is deliberately lane-agnostic.  It takes label counts or a label
vector and nothing else, so it applies to any endpoint in the benchmark that
scores a ranking, not only to the cohort it was first written for.  It uses
only the standard library and the package's own pure-Python metrics, so it runs
anywhere the package imports.
"""

from __future__ import annotations

from dataclasses import dataclass
from random import Random
from typing import Sequence

from .metrics import MetricError, average_precision, spearman_correlation


#: One-line rule this module exists to enforce.
GUIDANCE = (
    "Report a ranking metric against a random-score reference computed on the "
    "same label vector, never against class prevalence, and take the p-value "
    "from a permutation of the model's own scores so the model's ties are "
    "respected. Where the evaluation set has declared strata -- study, cohort, "
    "batch, platform, genomic block -- confine that permutation to within each "
    "stratum and report the stratified and global nulls side by side: a global "
    "shuffle destroys the nuisance structure along with the signal, crediting "
    "the observed statistic for structure the null does not have."
)

DEFAULT_DRAWS = 10_000
DEFAULT_SEED = 20260825


@dataclass(frozen=True)
class NullReference:
    """Summary of a metric's distribution under a stated null."""

    metric: str
    null: str
    mean: float
    sd: float
    percentile_95: float
    percentile_99: float
    maximum: float
    n_draws: int
    seed: int
    n_positive: int
    n_negative: int
    prevalence: float
    distinct_score_values: int | None
    percentile_scale: str = "signed_one_sided"
    #: Number of strata the permutation was confined to, or ``None`` when the
    #: null was drawn globally.  Both stratum fields are defaulted so existing
    #: construction sites and ``to_dict`` consumers are unaffected.
    n_strata: int | None = None
    #: Strata that cannot destroy any label-score association because they hold
    #: fewer than two members or only one distinct label.  Reported rather than
    #: dropped: a null built mostly from these is weaker than it looks.
    non_contributing_strata: int = 0

    @property
    def bias_over_prevalence(self) -> float:
        """How far a random scorer sits above the quoted prevalence baseline."""

        return self.mean - self.prevalence

    def to_dict(self) -> dict[str, object]:
        return {
            "metric": self.metric,
            "null": self.null,
            "percentile_scale": self.percentile_scale,
            "null_mean": self.mean,
            "null_sd": self.sd,
            "null_percentile_95": self.percentile_95,
            "null_percentile_99": self.percentile_99,
            "null_maximum": self.maximum,
            "n_draws": self.n_draws,
            "seed": self.seed,
            "n_positive": self.n_positive,
            "n_negative": self.n_negative,
            "prevalence": self.prevalence,
            "distinct_score_values": self.distinct_score_values,
            "bias_over_prevalence": self.bias_over_prevalence,
            "n_strata": self.n_strata,
            "non_contributing_strata": self.non_contributing_strata,
        }


def _percentile(ordered: Sequence[float], probability: float) -> float:
    if not ordered:
        raise MetricError("no draws to take a percentile from")
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _summarise(
    draws: list[float],
    *,
    metric: str,
    null: str,
    n_positive: int,
    n_negative: int,
    seed: int,
    distinct_score_values: int | None,
    absolute: bool = False,
    n_strata: int | None = None,
    non_contributing_strata: int = 0,
) -> NullReference:
    total = n_positive + n_negative
    mean = sum(draws) / len(draws)
    variance = sum((value - mean) ** 2 for value in draws) / (len(draws) - 1)
    ordered = sorted(abs(value) for value in draws) if absolute else sorted(draws)
    return NullReference(
        percentile_scale="absolute_two_sided" if absolute else "signed_one_sided",
        metric=metric,
        null=null,
        mean=mean,
        sd=variance**0.5,
        percentile_95=_percentile(ordered, 0.95),
        percentile_99=_percentile(ordered, 0.99),
        maximum=ordered[-1],
        n_draws=len(draws),
        seed=seed,
        n_positive=n_positive,
        n_negative=n_negative,
        prevalence=n_positive / total if total else float("nan"),
        distinct_score_values=distinct_score_values,
        n_strata=n_strata,
        non_contributing_strata=non_contributing_strata,
    )


def _scores(rng: Random, size: int, distinct_score_values: int | None) -> list[float]:
    if distinct_score_values is None:
        return [rng.random() for _ in range(size)]
    if distinct_score_values < 2:
        raise MetricError("a scorer needs at least two distinct values to rank")
    return [float(rng.randrange(distinct_score_values)) for _ in range(size)]


def random_score_reference(
    n_positive: int,
    n_negative: int,
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = DEFAULT_SEED,
    distinct_score_values: int | None = None,
) -> NullReference:
    """Average precision reached by a random scorer on this label vector.

    ``distinct_score_values`` models a scorer that emits only that many levels.
    Leave it ``None`` for the continuous case, which is the detection floor to
    quote beside a reported AUPRC.
    """

    if n_positive < 1 or n_negative < 1:
        raise MetricError("a ranking null needs both classes present")
    labels = [1] * n_positive + [0] * n_negative
    rng = Random(seed)
    draws: list[float] = []
    while len(draws) < n_draws:
        scores = _scores(rng, len(labels), distinct_score_values)
        draws.append(average_precision(labels, scores))
    return _summarise(
        draws,
        metric="average_precision",
        null=(
            "continuous_random_score"
            if distinct_score_values is None
            else f"random_score_with_{distinct_score_values}_distinct_values"
        ),
        n_positive=n_positive,
        n_negative=n_negative,
        seed=seed,
        distinct_score_values=distinct_score_values,
    )


def permutation_reference(
    labels: Sequence[int],
    scores: Sequence[float],
    *,
    n_permutations: int = DEFAULT_DRAWS,
    seed: int = DEFAULT_SEED,
) -> tuple[NullReference, float, float]:
    """Permute a model's own score vector, preserving its tie structure.

    Returns the null summary, the observed average precision, and the one-sided
    permutation p-value.  Use this for the p-value and
    ``random_score_reference`` for the detection floor; they answer different
    questions and neither substitutes for the other.
    """

    observed = average_precision(list(labels), list(scores))
    working = list(scores)
    rng = Random(seed)
    draws: list[float] = []
    for _ in range(n_permutations):
        rng.shuffle(working)
        draws.append(average_precision(list(labels), working))
    positives = sum(labels)
    reference = _summarise(
        draws,
        metric="average_precision",
        null="label_permutation_of_the_model_score_vector",
        n_positive=positives,
        n_negative=len(labels) - positives,
        seed=seed,
        distinct_score_values=len(set(scores)),
    )
    extreme = sum(1 for value in draws if value >= observed)
    return reference, observed, (extreme + 1) / (n_permutations + 1)


def _stratum_groups(
    strata: Sequence[object],
    length: int,
    name: str,
) -> dict[object, list[int]]:
    if len(strata) != length:
        raise MetricError(f"{name} must align with the observation vector")
    groups: dict[object, list[int]] = {}
    for index, stratum in enumerate(strata):
        groups.setdefault(stratum, []).append(index)
    if not groups:
        raise MetricError("no strata supplied")
    return groups


def stratified_permutation_reference(
    labels: Sequence[int],
    scores: Sequence[float],
    strata: Sequence[object],
    *,
    n_permutations: int = DEFAULT_DRAWS,
    seed: int = DEFAULT_SEED,
) -> tuple[NullReference, float, float]:
    """Permute the score vector *within* each stratum.

    A global permutation destroys the signal and the nuisance structure
    together, so an observed statistic is credited for between-stratum
    structure the null does not have.  Confining the shuffle to each stratum
    preserves that structure and destroys only the within-stratum association,
    which is the association a model is supposed to have learned.

    Report this beside :func:`permutation_reference` rather than instead of it;
    the gap between the two is itself the evidence that the evaluation set
    carries structure.

    Returns the null summary, the observed average precision, and the one-sided
    permutation p-value.
    """

    label_values = list(labels)
    score_values = list(scores)
    if len(label_values) != len(score_values):
        raise MetricError("labels and scores must align")
    groups = _stratum_groups(strata, len(label_values), "strata")

    permutable: list[list[int]] = []
    non_contributing = 0
    for indices in groups.values():
        distinct_labels = {label_values[index] for index in indices}
        if len(indices) < 2 or len(distinct_labels) < 2:
            # No label-score association exists inside this stratum for a
            # permutation to break.  Counted, never silently folded in.
            non_contributing += 1
        if len(indices) >= 2:
            permutable.append(indices)
    if non_contributing == len(groups):
        raise MetricError(
            "every stratum is non-contributing, so a within-stratum permutation "
            "destroys no association and the null collapses onto the observed "
            "value; widen the strata or use the global permutation reference"
        )

    observed = average_precision(label_values, score_values)
    rng = Random(seed)
    working = list(score_values)
    draws: list[float] = []
    for _ in range(n_permutations):
        for indices in permutable:
            block = [working[index] for index in indices]
            rng.shuffle(block)
            for index, value in zip(indices, block):
                working[index] = value
        draws.append(average_precision(label_values, working))
    positives = sum(label_values)
    reference = _summarise(
        draws,
        metric="average_precision",
        null="label_permutation_within_strata",
        n_positive=positives,
        n_negative=len(label_values) - positives,
        seed=seed,
        distinct_score_values=len(set(score_values)),
        n_strata=len(groups),
        non_contributing_strata=non_contributing,
    )
    extreme = sum(1 for value in draws if value >= observed)
    return reference, observed, (extreme + 1) / (n_permutations + 1)


def stratified_spearman_reference(
    outcome: Sequence[float],
    scores: Sequence[float],
    strata: Sequence[object],
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = DEFAULT_SEED,
) -> tuple[NullReference, float, float]:
    """Continuous analogue of :func:`stratified_permutation_reference`.

    Permutes the *outcome* within each stratum, so a correlation is judged
    against what a scorer reaches once between-stratum structure is held fixed.
    A stratum whose outcome is constant cannot contribute and is counted.
    """

    outcome_values = list(outcome)
    score_values = list(scores)
    if len(outcome_values) != len(score_values):
        raise MetricError("outcome and scores must align")
    if len(outcome_values) < 3:
        raise MetricError("a correlation null needs at least three observations")
    groups = _stratum_groups(strata, len(outcome_values), "strata")

    permutable: list[list[int]] = []
    non_contributing = 0
    for indices in groups.values():
        distinct_outcomes = {outcome_values[index] for index in indices}
        if len(indices) < 2 or len(distinct_outcomes) < 2:
            non_contributing += 1
        if len(indices) >= 2:
            permutable.append(indices)
    if non_contributing == len(groups):
        raise MetricError(
            "every stratum has a constant or single-member outcome, so a "
            "within-stratum permutation destroys no association; widen the "
            "strata or use the global reference"
        )

    observed = spearman_correlation(score_values, outcome_values)
    rng = Random(seed)
    working = list(outcome_values)
    draws: list[float] = []
    for _ in range(n_draws):
        for indices in permutable:
            block = [working[index] for index in indices]
            rng.shuffle(block)
            for index, value in zip(indices, block):
                working[index] = value
        draws.append(spearman_correlation(score_values, working))
    reference = _summarise(
        draws,
        metric="absolute_spearman",
        null="outcome_permutation_within_strata",
        n_positive=0,
        n_negative=len(outcome_values),
        seed=seed,
        distinct_score_values=len(set(score_values)),
        absolute=True,
        n_strata=len(groups),
        non_contributing_strata=non_contributing,
    )
    extreme = sum(1 for value in draws if abs(value) >= abs(observed))
    return reference, observed, (extreme + 1) / (n_draws + 1)


def spearman_reference(
    outcome: Sequence[float],
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = DEFAULT_SEED,
) -> NullReference:
    """Absolute Spearman correlation a continuous random scorer reaches.

    Ties in the *outcome* matter here in the same way ties in the scores matter
    for average precision, so the reference is computed against the real
    outcome vector rather than against a rank-uniform idealisation.
    """

    values = list(outcome)
    if len(values) < 3:
        raise MetricError("a correlation null needs at least three observations")
    rng = Random(seed)
    draws: list[float] = []
    while len(draws) < n_draws:
        scores = [rng.random() for _ in values]
        draws.append(spearman_correlation(scores, values))
    return _summarise(
        draws,
        metric="absolute_spearman",
        null="continuous_random_score",
        n_positive=0,
        n_negative=len(values),
        seed=seed,
        distinct_score_values=None,
        absolute=True,
    )


def bias_table(
    configurations: Sequence[tuple[int, int]],
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = DEFAULT_SEED,
) -> list[dict[str, object]]:
    """Tabulate the prevalence-versus-random-scorer gap across cohort shapes."""

    rows: list[dict[str, object]] = []
    for n_positive, n_negative in configurations:
        reference = random_score_reference(
            n_positive, n_negative, n_draws=n_draws, seed=seed
        )
        rows.append(
            {
                "n": n_positive + n_negative,
                "n_positive": n_positive,
                "n_negative": n_negative,
                "prevalence": reference.prevalence,
                "random_score_mean": reference.mean,
                "bias_over_prevalence": reference.bias_over_prevalence,
                "random_score_percentile_95": reference.percentile_95,
                "random_score_percentile_99": reference.percentile_99,
                "prevalence_understates_the_floor_by": (
                    reference.percentile_95 - reference.prevalence
                ),
            }
        )
    return rows


def tie_structure_table(
    n_positive: int,
    n_negative: int,
    distinct_value_grid: Sequence[int | None],
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = DEFAULT_SEED,
) -> list[dict[str, object]]:
    """Tabulate how a scorer's tie structure narrows its own null."""

    rows: list[dict[str, object]] = []
    for distinct in distinct_value_grid:
        reference = random_score_reference(
            n_positive,
            n_negative,
            n_draws=n_draws,
            seed=seed,
            distinct_score_values=distinct,
        )
        rows.append(
            {
                "distinct_score_values": (
                    "continuous" if distinct is None else distinct
                ),
                "null_mean": reference.mean,
                "null_sd": reference.sd,
                "null_percentile_95": reference.percentile_95,
                "null_percentile_99": reference.percentile_99,
            }
        )
    return rows
