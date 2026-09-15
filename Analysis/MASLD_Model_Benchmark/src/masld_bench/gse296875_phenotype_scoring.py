"""Pooled-only scoring for the GSE296875 donor-by-lineage phenotype campaign.

Fold 3 holds four donors and zero fibrosis positives, so a within-fold AUPRC is
undefined there.  A mean of per-fold metrics is also a different estimand than a
metric computed once on pooled out-of-fold predictions, whether or not every
fold happens to be defined.  This module makes the per-fold path unreachable
rather than discouraged:

1.  ``load_predictions`` projects the frozen prediction table to ``donor_id``
    and ``prediction``.  The fold column is never bound to a name, so no
    downstream caller has fold information to group by.
2.  ``registered_roster`` serves rosters only for names in a closed registry.
    No per-fold scope exists, and none can be constructed at call time.
3.  ``pooled_metric`` requires exact roster equality.  One fold's donors are a
    proper subset and raise; an unmasked donor is a superset and also raises.
4.  No per-fold metric function is defined here, so there is nothing to call by
    accident.

The campaign TaskSpec asks for a ``paired_donor_cluster_bootstrap``.  That
requirements are satisfied by ``pooled_donor_bootstrap``: it is paired and it
resamples donors.  What does not apply is the existing
``evaluators.stats.paired_cluster_bootstrap`` implementation, which reduces
within each cluster and then averages the differences.  That algorithm assumes
an additive metric, and neither pooled Spearman nor pooled average precision is
additive over donors.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from random import Random
from typing import Iterable, Mapping, Sequence

from .evaluators.metrics import average_precision, spearman_correlation


class PooledMetricViolation(RuntimeError):
    """Raised when a metric is asked for anything other than a full roster."""


class UnknownScope(RuntimeError):
    """Raised for any scope outside the closed registry."""


class PredictionTableError(RuntimeError):
    """Raised when the frozen prediction table does not match its requirements."""


ENDPOINTS = ("steatosis", "fibrosis")

METRIC_BY_ENDPOINT = {
    "steatosis": "spearman",
    "fibrosis": "auprc",
}

WELLS = ("well1", "well2", "well3", "well4", "well5", "well6", "well7", "well8")

#: The frozen fragment-membership primary partition.  Secondary lineages
#: (endothelial_cell, b_cell) are reported but are never family members.
PRIMARY_LINEAGES_DEFAULT = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)

#: Closed scope registry.  There is deliberately no per-fold scope here, and
#: ``registered_roster`` refuses any name absent from this set.
SCOPES = frozenset(
    ("pooled_all_donors", "adult_only_scored", "adult_only_refit")
    + tuple(f"leave_one_well_out_{well}" for well in WELLS)
)

#: Scopes that are confirmatory.  Everything else is a labelled sensitivity and
#: never enters the Benjamini-Hochberg family.
CONFIRMATORY_SCOPES = frozenset({"pooled_all_donors"})

PREDICTION_FIELDS = (
    "donor_id",
    "arm",
    "lineage_scope",
    "endpoint",
    "prediction",
    "outer_fold",
)


@dataclass(frozen=True)
class EndpointRoster:
    """The exact donor set a pooled metric must be computed over."""

    endpoint: str
    scope: str
    donor_ids: frozenset[str]

    def __post_init__(self) -> None:
        if self.endpoint not in ENDPOINTS:
            raise UnknownScope(f"unknown endpoint: {self.endpoint}")
        if self.scope not in SCOPES:
            raise UnknownScope(f"unregistered scope: {self.scope}")
        if len(self.donor_ids) < 2:
            raise PooledMetricViolation(
                f"{self.scope}/{self.endpoint} roster is degenerate"
            )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def registered_roster(
    endpoint: str,
    scope: str,
    endpoint_rows: Sequence[Mapping[str, str]],
    fold_rows: Sequence[Mapping[str, str]],
) -> EndpointRoster:
    """Build the donor roster for one registered scope from frozen tables.

    The roster comes from the frozen endpoint mask, never from the prediction
    table, so a model that silently failed to emit a donor cannot shrink the
    denominator it is scored against.
    """

    if endpoint not in ENDPOINTS:
        raise UnknownScope(f"unknown endpoint: {endpoint}")
    if scope not in SCOPES:
        raise UnknownScope(f"unregistered scope: {scope}")

    observed_field = f"{endpoint}_observed"
    observed = {
        row["donor_id"]
        for row in endpoint_rows
        if row[observed_field] == "true"
    }
    adults = {row["donor_id"] for row in endpoint_rows if row["is_adult"] == "true"}
    well_by_donor = {row["donor_id"]: row["well_id"] for row in fold_rows}

    if scope == "pooled_all_donors":
        donors = observed
    elif scope in {"adult_only_scored", "adult_only_refit"}:
        donors = observed & adults
    else:
        well = scope.removeprefix("leave_one_well_out_")
        donors = {donor for donor in observed if well_by_donor.get(donor) == well}

    return EndpointRoster(
        endpoint=endpoint, scope=scope, donor_ids=frozenset(donors)
    )


def load_predictions(
    path: Path,
    *,
    arm: str,
    lineage_scope: str,
    endpoint: str,
) -> dict[str, float]:
    """Project the frozen prediction table to ``donor_id -> prediction``.

    ``outer_fold`` is required to be present in the file, so the output file stays
    auditable, but it is never bound to a name.  The mapping returned here
    carries no fold information, which is what makes a per-fold split
    inexpressible downstream rather than merely discouraged.
    """

    rows = read_tsv(path)
    if not rows:
        raise PredictionTableError(f"prediction table is empty: {path}")
    missing = set(PREDICTION_FIELDS) - set(rows[0])
    if missing:
        raise PredictionTableError(
            f"prediction table is missing columns: {sorted(missing)}"
        )
    selected: dict[str, float] = {}
    for row in rows:
        if (
            row["arm"] != arm
            or row["lineage_scope"] != lineage_scope
            or row["endpoint"] != endpoint
        ):
            continue
        donor = row["donor_id"]
        if donor in selected:
            raise PredictionTableError(
                f"duplicate out-of-fold prediction for donor {donor}"
            )
        selected[donor] = float(row["prediction"])
    if not selected:
        raise PredictionTableError(
            f"no rows for {arm}/{lineage_scope}/{endpoint} in {path}"
        )
    return selected


def endpoint_values(
    endpoint: str, endpoint_rows: Sequence[Mapping[str, str]]
) -> dict[str, float]:
    """Read observed outcome values keyed by donor, honouring the frozen mask."""

    if endpoint not in ENDPOINTS:
        raise UnknownScope(f"unknown endpoint: {endpoint}")
    values: dict[str, float] = {}
    for row in endpoint_rows:
        if row[f"{endpoint}_observed"] != "true":
            continue
        if endpoint == "steatosis":
            values[row["donor_id"]] = float(row["steatosis_numeric"])
        else:
            values[row["donor_id"]] = 1.0 if row["fibrosis_any"] == "true" else 0.0
    return values


def _metric(kind: str, labels: Sequence[float], scores: Sequence[float]) -> float:
    if kind == "spearman":
        return spearman_correlation(scores, labels)
    if kind == "auprc":
        return average_precision([int(value) for value in labels], scores)
    raise ValueError(f"unknown metric kind: {kind}")


def _degenerate(kind: str, labels: Sequence[float], scores: Sequence[float]) -> bool:
    if kind == "auprc":
        positives = sum(1 for value in labels if value > 0.0)
        return positives == 0 or positives == len(labels)
    return len(set(labels)) < 2 or len(set(scores)) < 2


def pooled_metric(
    kind: str,
    predictions: Mapping[str, float],
    outcomes: Mapping[str, float],
    roster: EndpointRoster,
) -> float:
    """Compute one metric once, over exactly the roster's donors.

    Equality is required rather than containment.  A single fold's donors are a
    proper subset and raise here; a donor the frozen mask excludes is a superset
    and also raises.
    """

    supplied = frozenset(predictions)
    if supplied != roster.donor_ids:
        short = sorted(roster.donor_ids - supplied)
        extra = sorted(supplied - roster.donor_ids)
        raise PooledMetricViolation(
            f"{roster.scope}/{roster.endpoint} must be scored on exactly "
            f"{len(roster.donor_ids)} donors, received {len(supplied)}; "
            f"missing={short} unexpected={extra}"
        )
    missing_outcomes = sorted(roster.donor_ids - set(outcomes))
    if missing_outcomes:
        raise PooledMetricViolation(
            f"{roster.scope}/{roster.endpoint} outcome missing for "
            f"{missing_outcomes}"
        )
    donors = sorted(roster.donor_ids)
    labels = [outcomes[donor] for donor in donors]
    scores = [predictions[donor] for donor in donors]
    if _degenerate(kind, labels, scores):
        raise PooledMetricViolation(
            f"{roster.scope}/{roster.endpoint} is degenerate on the full roster"
        )
    return _metric(kind, labels, scores)


@dataclass(frozen=True)
class BootstrapSummary:
    estimate: float
    lower: float
    upper: float
    n_donors: int
    n_resamples: int
    degenerate_resamples: int
    seed: int

    def to_dict(self) -> dict[str, object]:
        return {
            "estimate": self.estimate,
            "ci_lower": self.lower,
            "ci_upper": self.upper,
            "n_donors": self.n_donors,
            "n_resamples": self.n_resamples,
            "degenerate_resamples": self.degenerate_resamples,
            "seed": self.seed,
        }


def _percentile(sorted_values: Sequence[float], probability: float) -> float:
    if not sorted_values:
        raise PooledMetricViolation("no non-degenerate bootstrap replicates")
    position = probability * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = position - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def pooled_donor_bootstrap(
    kind: str,
    predictions: Mapping[str, float],
    outcomes: Mapping[str, float],
    roster: EndpointRoster,
    *,
    baseline: Mapping[str, float] | None = None,
    n_resamples: int = 10_000,
    confidence_level: float = 0.95,
    seed: int = 20260825,
) -> BootstrapSummary:
    """Paired bootstrap clustered on donors, recomputing the pooled metric.

    Donors are the resampling unit, so the interval reflects n=39 and not the
    number of donor-by-lineage units.  When ``baseline`` is supplied the same
    resampled donor multiset scores both arms, which is what makes the
    difference paired.

    Degenerate replicates are counted and reported, never silently dropped.
    """

    point = pooled_metric(kind, predictions, outcomes, roster)
    if baseline is not None:
        point -= pooled_metric(kind, baseline, outcomes, roster)

    donors = sorted(roster.donor_ids)
    rng = Random(seed)
    estimates: list[float] = []
    degenerate = 0
    for _ in range(n_resamples):
        drawn = [donors[rng.randrange(len(donors))] for _ in donors]
        labels = [outcomes[donor] for donor in drawn]
        scores = [predictions[donor] for donor in drawn]
        if _degenerate(kind, labels, scores):
            degenerate += 1
            continue
        value = _metric(kind, labels, scores)
        if baseline is not None:
            base_scores = [baseline[donor] for donor in drawn]
            if _degenerate(kind, labels, base_scores):
                degenerate += 1
                continue
            value -= _metric(kind, labels, base_scores)
        estimates.append(value)

    estimates.sort()
    tail = (1.0 - confidence_level) / 2.0
    return BootstrapSummary(
        estimate=point,
        lower=_percentile(estimates, tail),
        upper=_percentile(estimates, 1.0 - tail),
        n_donors=len(donors),
        n_resamples=n_resamples,
        degenerate_resamples=degenerate,
        seed=seed,
    )


@dataclass(frozen=True)
class PermutationNull:
    """A null distribution together with the side it is measured on.

    ``percentile_scale`` is not decoration. Average precision is one sided, so
    its percentiles are taken on the signed statistic. A rank correlation is
    two sided, so its percentiles must be taken on the absolute statistic to
    match the two-sided p-value reported beside them. Two objects in one record
    once carried the same field name meaning different sides; this field is what
    makes that impossible to repeat.
    """

    mean: float
    sd: float
    percentile_95: float
    percentile_99: float
    maximum: float
    n_permutations: int
    seed: int
    observed: float
    p_value: float
    percentile_scale: str

    def to_dict(self) -> dict[str, object]:
        return {
            "percentile_scale": self.percentile_scale,
            "null_mean": self.mean,
            "null_sd": self.sd,
            "null_percentile_95": self.percentile_95,
            "null_percentile_99": self.percentile_99,
            "null_maximum": self.maximum,
            "n_permutations": self.n_permutations,
            "seed": self.seed,
            "observed": self.observed,
            "permutation_p_value": self.p_value,
        }


def permutation_null(
    kind: str,
    predictions: Mapping[str, float],
    outcomes: Mapping[str, float],
    roster: EndpointRoster,
    *,
    n_permutations: int = 10_000,
    seed: int = 20260825,
) -> PermutationNull:
    """Random-score reference for the metric on this exact roster.

    Average precision is upward biased at this sample size, so the prevalence
    line is not the right null for the fibrosis endpoint.  The comparison that
    belongs next to every reported AUPRC is this permutation distribution.
    """

    observed = pooled_metric(kind, predictions, outcomes, roster)
    donors = sorted(roster.donor_ids)
    labels = [outcomes[donor] for donor in donors]
    scores = [predictions[donor] for donor in donors]
    rng = Random(seed)
    draws: list[float] = []
    for _ in range(n_permutations):
        shuffled = list(scores)
        rng.shuffle(shuffled)
        draws.append(_metric(kind, labels, shuffled))
    mean = sum(draws) / len(draws)
    variance = sum((value - mean) ** 2 for value in draws) / (len(draws) - 1)
    # The percentile must be taken on the same side as the p-value below it.
    if kind == "spearman":
        ordered = sorted(abs(value) for value in draws)
        scale = "absolute_two_sided"
        extreme = sum(1 for value in draws if abs(value) >= abs(observed))
    else:
        ordered = sorted(draws)
        scale = "signed_one_sided"
        extreme = sum(1 for value in draws if value >= observed)
    return PermutationNull(
        mean=mean,
        sd=variance**0.5,
        percentile_95=_percentile(ordered, 0.95),
        percentile_99=_percentile(ordered, 0.99),
        maximum=ordered[-1],
        n_permutations=n_permutations,
        seed=seed,
        observed=observed,
        p_value=(extreme + 1) / (n_permutations + 1),
        percentile_scale=scale,
    )


def random_score_reference(
    kind: str,
    outcomes: Mapping[str, float],
    roster: EndpointRoster,
    *,
    n_draws: int = 10_000,
    seed: int = 20260825,
) -> PermutationNull:
    """Reference distribution for a *continuous* random scorer on this roster.

    This is the null that belongs in the power statement, because all three
    arms emit continuous scores.  It differs from :func:`permutation_null`,
    which permutes a specific model's own prediction vector and therefore
    inherits that vector's tie structure: a heavily tied scorer has a visibly
    tighter null than a continuous one.  Report both and do not substitute one
    for the other.

    ``observed`` and ``p_value`` are not meaningful without a model, so they are
    reported as NaN.
    """

    donors = sorted(roster.donor_ids)
    labels = [outcomes[donor] for donor in donors]
    rng = Random(seed)
    draws: list[float] = []
    while len(draws) < n_draws:
        scores = [rng.random() for _ in donors]
        if _degenerate(kind, labels, scores):
            continue
        draws.append(_metric(kind, labels, scores))
    mean = sum(draws) / len(draws)
    variance = sum((value - mean) ** 2 for value in draws) / (len(draws) - 1)
    if kind == "spearman":
        ordered = sorted(abs(value) for value in draws)
        scale = "absolute_two_sided"
    else:
        ordered = sorted(draws)
        scale = "signed_one_sided"
    return PermutationNull(
        mean=mean,
        sd=variance**0.5,
        percentile_95=_percentile(ordered, 0.95),
        percentile_99=_percentile(ordered, 0.99),
        maximum=ordered[-1],
        n_permutations=n_draws,
        seed=seed,
        observed=float("nan"),
        p_value=float("nan"),
        percentile_scale=scale,
    )


def confirmatory_family(
    lineages: Iterable[str],
) -> tuple[tuple[str, str, str], ...]:
    """Enumerate the Benjamini-Hochberg family: molecular arm, pooled scope.

    Two endpoints times the ``all_lineage`` scope plus each frozen primary
    lineage.  Metadata-only and combined arms are comparators reported as
    paired differences, not independent hypotheses, so they are absent here.
    """

    scopes = ("all_lineage", *sorted(lineages))
    return tuple(
        ("molecular", scope, endpoint)
        for endpoint in ENDPOINTS
        for scope in scopes
    )
