"""Answer the frozen axis-count question, with the prespecification proved first.

Stage 0b of the MASLD showcase model, job 2 of 2. The wrapper verifies the
frozen prespecification tree with ``sha256sum --check --strict`` before this
runs; this script then re-reads the criteria and their literal thresholds *out
of that frozen file* rather than restating them, so a threshold cannot drift
between the freeze and the answer without aborting the run.

The question is how many histological outcome axes this substrate supports.
Stage 0 answered a different question -- whether the three NAS components are
three separable axes -- and returned STOP. Stage 0b was derived after seeing
that answer; the prespecification says so in those words and cites Stage 0 by
digest.

D1, decisive. Per-gene partial Spearman with the NAS sum adjusting for fibrosis
stage, and with fibrosis stage adjusting for the NAS sum. A direction survives
only if its BH 0.05 count is above zero *and* above the 95th percentile of a
matched permutation count null built from the same residuals. Both directions
must survive for there to be more than one axis.

D2 decides two axes against three: steatosis adjusting for
activity{ballooning + lobular inflammation}, and the reverse.

D3 is a diagnostic and never a check: genes reaching BH for a composite while
reaching BH for none of its constituents alone. It is confounded by power in a
direction that always favours the composite, so it is reported with the tie
structure of every axis beside it and with the reverse count.

Nulls are measured, never assumed. The count null permutes the
covariate-residualized exposure and re-residualizes it, then recomputes the
whole family, so every draw yields a complete BH count under the same tie
structure the observed count was computed with. The primary reference permutes
within ``recorded_sex``, which Stage 0 measured to be a genuine stratum here;
the global permutation is reported beside it. ``outer_fold`` is not a stratum
and Stage 0 ruled it out correctly.

``lobular_necrosis`` is deposited but is not a NASH-CRN NAS component. It enters
no composite and no criterion.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.stats import rankdata
from scipy.special import stdtr
from scipy.stats import t as student_t

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.evaluators.auprc_reference import (
    spearman_reference,
    stratified_spearman_reference,
)
from masld_bench.evaluators.metrics import _average_ranks, spearman_correlation
from masld_bench.evaluators.stats import benjamini_hochberg, gain_concentration


class AxisCountError(RuntimeError):
    """Raised when the substrate cannot support the frozen question."""


COMPONENTS = ("steatosis", "ballooning", "lobular_inflammation")
NAS_ACTIVITY = ("steatosis", "ballooning", "lobular_inflammation")
SAF_ACTIVITY = ("ballooning", "lobular_inflammation")

D1_FAMILIES = (
    ("nas_activity_sum", "fibrosis"),
    ("fibrosis", "nas_activity_sum"),
)
D2_FAMILIES = (
    ("steatosis", "saf_activity_sum"),
    ("saf_activity_sum", "steatosis"),
)

BH_ALPHA = 0.05
COUNT_FLOOR_PERCENTILE = 95.0

FROZEN_THRESHOLDS = {
    "d1_activity_and_fibrosis_are_independent": (
        "in each direction, BH 0.05 count > 0 and above the 95th percentile of "
        "the matched within-sex permutation count null"
    ),
    "d2_steatosis_and_saf_activity_are_independent": (
        "in each direction, BH 0.05 count > 0 and above the 95th percentile of "
        "the matched within-sex permutation count null"
    ),
}


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def read_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def load_prespecification(directory: Path) -> dict:
    """Read the frozen criteria and refuse to run against a drifted freeze."""

    verify_frozen_tree(directory)
    payload = json.loads(
        (directory / "axis_count_prespecification.json").read_text(encoding="utf-8")
    )
    criteria = payload["criteria"]
    if sorted(criteria) != sorted(FROZEN_THRESHOLDS):
        raise AxisCountError("the frozen criteria are not the two expected ones")
    for name, expected in FROZEN_THRESHOLDS.items():
        if criteria[name]["threshold"] != expected:
            raise AxisCountError(
                f"{name} threshold drifted since the freeze: "
                f"{criteria[name]['threshold']!r} is not {expected!r}"
            )
    if criteria["d1_activity_and_fibrosis_are_independent"].get("decisive") is not True:
        raise AxisCountError("d1 lost its decisive mark since the freeze")
    provenance = payload["honest_provenance"]
    for key in (
        "this_question_was_derived_after_seeing_stage_0s_answer",
        "criteria_fixed_before_any_expression_value_was_read",
    ):
        if provenance.get(key) is not True:
            raise AxisCountError(f"the frozen provenance does not assert {key}")
    return payload


# --------------------------------------------------------------------------
# ranks and residuals
# --------------------------------------------------------------------------


def average_ranks(matrix: np.ndarray) -> np.ndarray:
    """Average ranks down each column, matching the package's tie handling."""

    return rankdata(matrix, method="average", axis=0)


def validate_rank_agreement(
    matrix: np.ndarray, columns: Sequence[int]
) -> dict[str, object]:
    """Prove the vectorized ranks equal the package's own, on real columns."""

    ranked = average_ranks(matrix[:, list(columns)])
    worst = 0.0
    for offset, column in enumerate(columns):
        reference = np.asarray(_average_ranks(list(matrix[:, column])), dtype=float)
        worst = max(worst, float(np.max(np.abs(ranked[:, offset] - reference))))
    return {
        "columns_checked": len(columns),
        "max_absolute_rank_difference": worst,
        "agrees_with_masld_bench_average_ranks": worst == 0.0,
    }


def centred_ranks(values: np.ndarray) -> np.ndarray:
    ranked = rankdata(np.asarray(values, dtype=float), method="average").astype(float)
    return ranked - ranked.mean()


def unit_ranks(values: np.ndarray) -> np.ndarray:
    centred = centred_ranks(values)
    norm = float(np.sqrt((centred**2).sum()))
    if norm == 0.0:
        raise AxisCountError("a constant axis reached the correlation step")
    return centred / norm


def residualize(target: np.ndarray, unit_covariate: np.ndarray) -> np.ndarray:
    """Remove the covariate's linear-on-ranks component, column-wise.

    ``unit_covariate`` is already centred and unit-norm, so the projection is a
    single outer product and the intercept is handled by the centring.
    """

    if target.ndim == 1:
        return target - unit_covariate * float(unit_covariate @ target)
    return target - np.outer(unit_covariate, unit_covariate @ target)


def unit_columns(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Scale each column to unit norm, returning the norms so zeros are visible."""

    norms = np.sqrt((matrix**2).sum(axis=0))
    keep = norms > 0.0
    scaled = np.zeros_like(matrix)
    scaled[:, keep] = matrix[:, keep] / norms[keep]
    return scaled, norms


# --------------------------------------------------------------------------
# per-gene inference
# --------------------------------------------------------------------------


def t_p_values(correlations: np.ndarray, residual_df: int) -> np.ndarray:
    """Two-sided Student t p-values for a partial correlation, df = n - 3.

    Reported beside a permutation p from the same draws. The count null exists
    precisely because this asymptotic form is not guaranteed calibrated at this
    tie structure.
    """

    bounded = np.clip(np.abs(correlations), 0.0, 1.0 - 1e-12)
    statistic = bounded * math.sqrt(residual_df) / np.sqrt(1.0 - bounded**2)
    return 2.0 * (1.0 - stdtr(residual_df, statistic))


def bh_count(p_values: np.ndarray, alpha: float = BH_ALPHA) -> int:
    """Number of hypotheses Benjamini-Hochberg rejects at ``alpha``.

    The step-up rule written out directly rather than by adjusting every
    p-value, because the null draws need this many thousands of times. Proved
    against ``masld_bench.evaluators.stats.benjamini_hochberg`` on the observed
    vector before any draw is taken.
    """

    m = p_values.size
    if m == 0:
        return 0
    ordered = np.sort(p_values)
    thresholds = alpha * np.arange(1, m + 1) / m
    passing = np.flatnonzero(ordered <= thresholds)
    return int(passing[-1] + 1) if passing.size else 0


def bh_critical_abs_r(m: int, residual_df: int, alpha: float = BH_ALPHA) -> np.ndarray:
    """The |partial r| a gene at BH rank k must reach, for k = 1..m.

    The two-sided t p-value is strictly decreasing in ``|r|`` at fixed df, so
    BH on those p-values is exactly a comparison of the descending-sorted
    ``|r|`` against this vector. Precomputing it once per family keeps the
    inverse-t out of the inner permutation loop, where it would otherwise be
    evaluated hundreds of millions of times.
    """

    ranks = np.arange(1, m + 1, dtype=float)
    critical_t = student_t.isf(alpha * ranks / (2.0 * m), residual_df)
    return critical_t / np.sqrt(critical_t**2 + residual_df)


def bh_count_from_abs_r(absolute: np.ndarray, critical: np.ndarray) -> int:
    """BH count read straight off the descending-sorted absolute statistic."""

    if absolute.size == 0:
        return 0
    ordered = np.sort(absolute)[::-1]
    passing = np.flatnonzero(ordered >= critical)
    return int(passing[-1] + 1) if passing.size else 0


def bh_counts_from_abs_r(block: np.ndarray, critical: np.ndarray) -> np.ndarray:
    """Vectorized :func:`bh_count_from_abs_r` over the rows of a null block."""

    ordered = -np.sort(-np.abs(block), axis=1)
    passing = ordered >= critical
    reversed_index = passing.shape[1] - 1 - np.argmax(passing[:, ::-1], axis=1)
    return np.where(passing.any(axis=1), reversed_index + 1, 0).astype(np.int64)


def critical_correlation(count: int, m: int, residual_df: int) -> float | None:
    """The smallest |partial r| a BH rejection at this count corresponds to."""

    if count <= 0 or m == 0:
        return None
    p = BH_ALPHA * count / m
    # invert the two-sided t p-value on |r|
    low, high = 0.0, 1.0 - 1e-12
    for _ in range(200):
        mid = 0.5 * (low + high)
        if float(t_p_values(np.asarray([mid]), residual_df)[0]) > p:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def _permuted_unit_exposure(
    residual_exposure: np.ndarray,
    unit_covariate: np.ndarray,
    order: np.ndarray,
) -> np.ndarray | None:
    """Permute the exposure residual, re-residualize, re-normalise.

    Re-residualizing matters: a raw permutation of the residual is no longer
    orthogonal to the covariate, so a gene that tracks the covariate would pick
    up spurious correlation. Returns ``None`` for the measure-zero draw that
    lands exactly on the covariate.
    """

    permuted = residualize(residual_exposure[order], unit_covariate)
    norm = float(np.sqrt((permuted**2).sum()))
    if norm == 0.0:
        return None
    return permuted / norm


def _stratum_orders(
    groups: Sequence[np.ndarray], generator: np.random.Generator, n: int
) -> np.ndarray:
    """One within-stratum permutation of the participant index."""

    order = np.arange(n)
    for indices in groups:
        if indices.size > 1:
            order[indices] = indices[generator.permutation(indices.size)]
    return order


def _stratum_index_groups(strata: Sequence[object]) -> dict[object, np.ndarray]:
    groups: dict[object, list[int]] = {}
    for index, value in enumerate(strata):
        groups.setdefault(value, []).append(index)
    return {key: np.asarray(value, dtype=int) for key, value in groups.items()}


def partial_family(
    unit_genes_residual: np.ndarray,
    residual_exposure: np.ndarray,
    unit_covariate: np.ndarray,
    *,
    exposure: str,
    covariate: str,
    n_participants: int,
    n_draws: int,
    block: int,
    seed: int,
    strata: Sequence[object] | None,
    null_name: str,
) -> dict[str, object]:
    """One direction: observed count, its BH family, and a matched count null.

    ``unit_genes_residual`` holds each gene's covariate-residualized rank
    vector scaled to unit norm, so a dot product with a unit exposure residual
    is exactly the rank-linear partial correlation.
    """

    residual_df = n_participants - 3
    norm = float(np.sqrt((residual_exposure**2).sum()))
    if norm == 0.0:
        raise AxisCountError(
            f"{exposure} is collinear with {covariate}; the partial is undefined"
        )
    unit_exposure = residual_exposure / norm

    observed = unit_genes_residual.T @ unit_exposure
    m = int(observed.size)
    critical = bh_critical_abs_r(m, residual_df)
    observed_p = t_p_values(observed, residual_df)
    count = bh_count(observed_p)
    if count != bh_count_from_abs_r(np.abs(observed), critical):
        raise AxisCountError("the fast BH path disagrees with the p-value path")

    generator = np.random.default_rng(seed)
    groups = (
        list(_stratum_index_groups(strata).values()) if strata is not None else None
    )
    draws = np.empty(n_draws, dtype=np.int64)
    maxima = np.empty(n_draws, dtype=float)
    exceedances = np.zeros(m, dtype=np.int64)
    absolute_observed = np.abs(observed)
    degenerate = 0
    drawn = 0
    while drawn < n_draws:
        size = min(block, n_draws - drawn)
        permuted = np.empty((size, n_participants), dtype=float)
        for index in range(size):
            order = (
                _stratum_orders(groups, generator, n_participants)
                if groups is not None
                else generator.permutation(n_participants)
            )
            candidate = _permuted_unit_exposure(
                residual_exposure, unit_covariate, order
            )
            if candidate is None:
                degenerate += 1
                candidate = unit_exposure
            permuted[index] = candidate
        absolute_block = np.abs(permuted @ unit_genes_residual)
        maxima[drawn : drawn + size] = absolute_block.max(axis=1)
        exceedances += (absolute_block >= absolute_observed).sum(axis=0)
        draws[drawn : drawn + size] = bh_counts_from_abs_r(absolute_block, critical)
        drawn += size

    permutation_p = (exceedances + 1) / (n_draws + 1)
    permutation_bh_count = bh_count(np.asarray(permutation_p, dtype=float))
    floor = float(np.percentile(draws, COUNT_FLOOR_PERCENTILE))
    extreme = int(np.sum(draws >= count))
    return {
        "exposure": exposure,
        "adjusted_for": covariate,
        "family": f"every gene in the realised {exposure} | {covariate} partial universe",
        "family_size": m,
        "residual_df": residual_df,
        "bh_alpha": BH_ALPHA,
        "genes_bh_below_0_05": int(count),
        "smallest_abs_partial_r_reaching_bh_0_05": (
            float(absolute_observed[observed_p <= BH_ALPHA * count / m].min())
            if count > 0
            else None
        ),
        "bh_critical_abs_partial_r_at_this_count": critical_correlation(
            count, m, residual_df
        ),
        "observed_max_abs_partial_r": float(absolute_observed.max()),
        "count_null": {
            "null": null_name,
            "n_draws": n_draws,
            "seed": seed,
            "degenerate_draws_reusing_the_observed_residual": degenerate,
            "null_mean_count": float(draws.mean()),
            "null_median_count": float(np.median(draws)),
            "null_percentile_95_count": floor,
            "null_percentile_99_count": float(np.percentile(draws, 99)),
            "null_max_count": int(draws.max()),
            "count_exceedance_p_value": (extreme + 1) / (n_draws + 1),
            "p_value_resolution_floor": 1.0 / (n_draws + 1),
            "observed_count_above_the_95th_percentile": bool(count > floor),
        },
        "familywise_floor": {
            "null": "residual_permutation_max_statistic_over_the_realised_family",
            "detectable_abs_partial_r_p95": float(np.percentile(maxima, 95)),
            "detectable_abs_partial_r_p99": float(np.percentile(maxima, 99)),
            "max_statistic_null_mean": float(maxima.mean()),
            "genes_above_the_familywise_floor": int(
                (absolute_observed > np.percentile(maxima, 95)).sum()
            ),
        },
        "permutation_p_sensitivity": {
            "why": (
                "the same draws give a per-gene permutation p, which does not "
                "assume the t form is calibrated at this tie structure"
            ),
            "genes_bh_below_0_05_on_permutation_p": int(permutation_bh_count),
            "p_value_resolution_floor": 1.0 / (n_draws + 1),
            "genes_at_the_permutation_p_floor": int(
                np.sum(np.asarray(permutation_p) <= 1.0 / (n_draws + 1))
            ),
            "bh_is_reachable_on_permutation_p_only_if_at_least_k_genes_sit_there": (
                math.ceil(m / (BH_ALPHA * (n_draws + 1)))
            ),
        },
    }


def marginal_family(
    unit_genes: np.ndarray,
    outcome: np.ndarray,
    *,
    axis: str,
    n_participants: int,
) -> dict[str, object]:
    """Marginal BH count for one axis, in the full non-constant universe."""

    observed = unit_genes.T @ unit_ranks(outcome)
    residual_df = n_participants - 2
    p_values = t_p_values(observed, residual_df)
    count = bh_count(p_values)
    return {
        "axis": axis,
        "family_size": int(observed.size),
        "residual_df": residual_df,
        "genes_bh_below_0_05": int(count),
        "observed_max_abs_spearman": float(np.abs(observed).max()),
        "bh_critical_abs_spearman_at_this_count": critical_correlation(
            count, int(observed.size), residual_df
        ),
    }


def tie_profile(values: np.ndarray, name: str) -> dict[str, object]:
    counts = Counter(int(value) for value in values.tolist())
    n = values.size
    tied = sum(count * (count - 1) // 2 for count in counts.values())
    return {
        "axis": name,
        "distinct_values": len(counts),
        "value_counts": {str(key): counts[key] for key in sorted(counts)},
        "largest_tied_block": max(counts.values()),
        "fraction_of_pairs_tied": tied / (n * (n - 1) // 2),
    }


# --------------------------------------------------------------------------
# vacuity
# --------------------------------------------------------------------------


def evaluate_gate(
    name: str, conditions: Sequence[dict[str, object]]
) -> dict[str, object]:
    """met over applicable, never met over total, and never a silent all([]).

    A condition is applicable only when its family could carry the statistic at
    all. A check with no applicable condition returns NO_APPLICABLE_CONDITIONS,
    which is not a pass and not a fail.
    """

    applicable = [entry for entry in conditions if entry["applicable"]]
    met = [entry for entry in applicable if entry["met"]]
    if not applicable:
        return {
            "gate": name,
            "verdict": "NO_APPLICABLE_CONDITIONS",
            "passed": False,
            "n_conditions_total": len(conditions),
            "n_conditions_applicable": 0,
            "n_conditions_met": 0,
            "reported_as": "met over applicable, never met over total",
            "why": (
                "every condition was inapplicable, so there is nothing to pass. "
                "all([]) is True and is guarded against here."
            ),
            "conditions": list(conditions),
        }
    return {
        "gate": name,
        "verdict": "MET" if len(met) == len(applicable) else "NOT_MET",
        "passed": len(met) == len(applicable),
        "n_conditions_total": len(conditions),
        "n_conditions_applicable": len(applicable),
        "n_conditions_met": len(met),
        "reported_as": "met over applicable, never met over total",
        "conditions": list(conditions),
    }


def direction_condition(result: dict[str, object]) -> dict[str, object]:
    count = int(result["genes_bh_below_0_05"])
    floor = float(result["count_null"]["null_percentile_95_count"])
    applicable = int(result["family_size"]) > 0
    met = bool(applicable and count > 0 and count > floor)
    return {
        "condition": f"{result['exposure']} | {result['adjusted_for']}",
        "applicable": applicable,
        "why_not_applicable": (
            None if applicable else "the realised partial family is empty"
        ),
        "family_size": int(result["family_size"]),
        "genes_bh_below_0_05": count,
        "count_null_percentile_95": floor,
        "count_exceedance_p_value": result["count_null"]["count_exceedance_p_value"],
        "met": met,
    }


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prespecification", type=Path, required=True)
    parser.add_argument("--endpoints", type=Path, required=True)
    parser.add_argument("--molecular", type=Path, required=True)
    parser.add_argument("--folds", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count-null-draws", type=int, default=20000)
    parser.add_argument("--null-block", type=int, default=250)
    parser.add_argument("--label-draws", type=int, default=20000)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260827)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise AxisCountError("refusing to overwrite an axis-count result")

    prespec = load_prespecification(arguments.prespecification)
    print("verified the frozen prespecification", flush=True)

    rows = read_table(arguments.endpoints)
    participants = [row["participant_id"] for row in rows]
    n = len(rows)
    if len(set(participants)) != n:
        raise AxisCountError("participant_id is not unique in the endpoint table")
    if n != prespec["cohort"]["n_participants"]:
        raise AxisCountError("the participant count differs from the freeze")

    folds = {
        row["participant_id"]: row["outer_fold"] for row in read_table(arguments.folds)
    }
    if folds != {row["participant_id"]: row["outer_fold"] for row in rows}:
        raise AxisCountError("the fold table disagrees with the endpoint table")

    components = {
        name: np.asarray([int(row[name]) for row in rows], dtype=float)
        for name in COMPONENTS
    }
    fibrosis = np.asarray([int(row["fibrosis"]) for row in rows], dtype=float)
    deposited_sum = np.asarray(
        [int(row["nash_crn_component_sum"]) for row in rows], dtype=float
    )
    nas_activity = sum(components[name] for name in NAS_ACTIVITY)
    saf_activity = sum(components[name] for name in SAF_ACTIVITY)
    if not np.array_equal(deposited_sum, nas_activity):
        raise AxisCountError("the deposited component sum is not the three components")
    necrosis_folded = np.array_equal(
        deposited_sum,
        nas_activity + np.asarray([int(row["lobular_necrosis"]) for row in rows], float),
    )
    if necrosis_folded:
        raise AxisCountError("lobular_necrosis cannot be inside the component sum")

    axes = {
        **components,
        "fibrosis": fibrosis,
        "nas_activity_sum": nas_activity,
        "saf_activity_sum": saf_activity,
    }
    ties = {name: tie_profile(values, name) for name, values in axes.items()}
    sexes = [row["recorded_sex"] for row in rows]

    # ---------------- label-level context, stratified reference primary ----
    criterion_axes = ("steatosis", "saf_activity_sum", "nas_activity_sum", "fibrosis")
    label_context = []
    for index, left in enumerate(criterion_axes):
        for right in criterion_axes[index + 1 :]:
            reference, observed, stratified_p = stratified_spearman_reference(
                axes[right].tolist(),
                axes[left].tolist(),
                sexes,
                n_draws=arguments.label_draws,
                seed=arguments.seed,
            )
            label_context.append({
                "pair": [left, right],
                "spearman": observed,
                "sex_stratified_null_primary": {
                    **reference.to_dict(),
                    "p_value_two_sided": stratified_p,
                    "p_value_resolution_floor": 1.0 / (arguments.label_draws + 1),
                },
                "continuous_random_scorer_reference": spearman_reference(
                    axes[right].tolist(),
                    n_draws=arguments.label_draws,
                    seed=arguments.seed,
                ).to_dict(),
                "inside_the_criteria": False,
            })
    print("label-level context done", flush=True)

    # ---------------- expression ----------------
    axis_rows = read_table(arguments.molecular / "participant_axis.tsv")
    if [row["participant_id"] for row in axis_rows] != participants:
        raise AxisCountError(
            "the molecular participant axis is not in the endpoint table's order"
        )
    features = read_table(arguments.molecular / "rna_feature_axis.tsv")
    gene_ids = [row["stable_gene_id"] for row in features]
    if len(set(gene_ids)) != len(gene_ids):
        raise AxisCountError("stable_gene_id is not unique on the feature axis")
    values = np.load(arguments.molecular / "rna_values.npy")
    if values.shape != (n, len(gene_ids)):
        raise AxisCountError("the expression matrix does not match its axes")
    print(f"loaded expression {values.shape}", flush=True)

    # The realised universe is re-derived here, never copied from Stage 0.
    ordered = np.sort(values, axis=0)
    realised = ((ordered[1:] != ordered[:-1]).sum(axis=0) + 1) >= 2
    del ordered
    universe = np.flatnonzero(realised)
    if universe.size == 0:
        raise AxisCountError("no gene varies across participants")
    realised_values = values[:, universe]
    print(
        f"realised gene universe: {universe.size} of {len(gene_ids)} "
        f"({len(gene_ids) - universe.size} constant across all {n} participants)",
        flush=True,
    )

    generator = np.random.default_rng(arguments.seed)
    sample = sorted(
        generator.choice(universe.size, size=min(300, universe.size), replace=False)
        .tolist()
    )
    rank_agreement = validate_rank_agreement(realised_values, sample)
    if not rank_agreement["agrees_with_masld_bench_average_ranks"]:
        raise AxisCountError("the vectorized ranks disagree with the package")

    gene_ranks = average_ranks(realised_values)
    gene_ranks -= gene_ranks.mean(axis=0, keepdims=True)
    unit_genes, marginal_norms = unit_columns(gene_ranks)
    if np.any(marginal_norms == 0.0):
        raise AxisCountError("a constant column survived the universe filter")

    # Prove the marginal association equals the package's own Spearman.
    marginal_check = []
    for name in ("steatosis", "fibrosis"):
        vector = unit_genes.T @ unit_ranks(axes[name])
        worst = 0.0
        for column in sample[:100]:
            direct = spearman_correlation(
                realised_values[:, column].tolist(), axes[name].tolist()
            )
            worst = max(worst, abs(direct - float(vector[column])))
        marginal_check.append({
            "axis": name,
            "genes_checked": 100,
            "max_absolute_difference_vs_masld_bench_spearman": worst,
        })
    if max(
        entry["max_absolute_difference_vs_masld_bench_spearman"]
        for entry in marginal_check
    ) > 1e-10:
        raise AxisCountError("the vectorized association disagrees with the package")

    # Prove the fast BH step-up equals the package's adjustment, once.
    probe = t_p_values(unit_genes.T @ unit_ranks(axes["steatosis"]), n - 2)
    package_count = int(
        np.sum(np.asarray(benjamini_hochberg(probe.tolist()), dtype=float) < BH_ALPHA)
    )
    if package_count != bh_count(probe):
        raise AxisCountError("the BH step-up disagrees with masld_bench")
    bh_agreement = {
        "checked_on": "the marginal steatosis family",
        "package_count": package_count,
        "step_up_count": bh_count(probe),
        "agrees_with_masld_bench_benjamini_hochberg": True,
    }
    print("rank, association and BH agreement validated", flush=True)

    # ---------------- the four partial families ----------------
    def build_family(
        exposure: str, covariate: str, strata: Sequence[object] | None, null_name: str,
        seed_offset: int,
    ) -> tuple[dict[str, object], dict[str, int]]:
        unit_covariate = unit_ranks(axes[covariate])
        residual_genes = residualize(gene_ranks, unit_covariate)
        scaled, norms = unit_columns(residual_genes)
        keep = norms > 0.0
        dropped = int((~keep).sum())
        result = partial_family(
            scaled[:, keep],
            residualize(centred_ranks(axes[exposure]), unit_covariate),
            unit_covariate,
            exposure=exposure,
            covariate=covariate,
            n_participants=n,
            n_draws=arguments.count_null_draws,
            block=arguments.null_block,
            seed=arguments.seed + seed_offset,
            strata=strata,
            null_name=null_name,
        )
        result["realised_universe"] = {
            "features_on_axis": len(gene_ids),
            "non_constant_across_all_participants": int(universe.size),
            "dropped_as_collinear_with_the_covariate": dropped,
            "entered_this_family": int(keep.sum()),
            "re_derived_not_copied": True,
        }
        return result, {"dropped": dropped}

    primary = {}
    global_null = {}
    for offset, (exposure, covariate) in enumerate(D1_FAMILIES + D2_FAMILIES):
        key = f"{exposure}|{covariate}"
        primary[key], _ = build_family(
            exposure, covariate, sexes,
            "within_recorded_sex_residual_permutation", 10 + offset,
        )
        print(f"primary family done: {key}", flush=True)
        global_null[key], _ = build_family(
            exposure, covariate, None,
            "global_residual_permutation", 50 + offset,
        )
        print(f"global-null family done: {key}", flush=True)

    # ---------------- checks ----------------
    d1 = evaluate_gate(
        "d1_activity_and_fibrosis_are_independent",
        [direction_condition(primary[f"{a}|{b}"]) for a, b in D1_FAMILIES],
    )
    d2 = evaluate_gate(
        "d2_steatosis_and_saf_activity_are_independent",
        [direction_condition(primary[f"{a}|{b}"]) for a, b in D2_FAMILIES],
    )

    d1_survivors = [
        entry["condition"] for entry in d1["conditions"]
        if entry["applicable"] and entry["met"]
    ]
    if d1["verdict"] == "NO_APPLICABLE_CONDITIONS" or (
        d2["verdict"] == "NO_APPLICABLE_CONDITIONS"
    ):
        outcome = "INDETERMINATE"
        grouping = "no applicable condition; the substrate cannot answer the question"
    elif d1["passed"] and d2["passed"]:
        outcome = "THREE_AXIS_SAF"
        grouping = "steatosis / activity{ballooning + lobular_inflammation} / fibrosis"
    elif d1["passed"]:
        outcome = "TWO_AXIS_ACTIVITY_FIBROSIS"
        grouping = "activity{full NAS sum} / fibrosis"
    elif d2["passed"]:
        outcome = "TWO_AXIS_STEATOSIS_ACTIVITY"
        grouping = (
            "steatosis / activity{ballooning + lobular_inflammation}, with "
            "fibrosis not independent of the NAS sum on this substrate"
        )
    elif d1_survivors == ["nas_activity_sum | fibrosis"]:
        outcome = "ONE_AXIS_ACTIVITY"
        grouping = "activity{full NAS sum} only"
    elif d1_survivors == ["fibrosis | nas_activity_sum"]:
        outcome = "ONE_AXIS_FIBROSIS"
        grouping = "fibrosis only"
    else:
        outcome = "ONE_AXIS_UNDETERMINED"
        grouping = "one axis; neither D1 direction survives, so neither names it"
    print(f"outcome: {outcome}", flush=True)

    # ---------------- diagnostics ----------------
    marginal = {
        name: marginal_family(unit_genes, values_vector, axis=name, n_participants=n)
        for name, values_vector in axes.items()
    }
    print("marginal counts done", flush=True)

    def bh_mask(outcome_vector: np.ndarray) -> np.ndarray:
        vector = unit_genes.T @ unit_ranks(outcome_vector)
        p_values = t_p_values(vector, n - 2)
        count = bh_count(p_values)
        if count == 0:
            return np.zeros(vector.size, dtype=bool)
        return p_values <= BH_ALPHA * count / p_values.size

    masks = {name: bh_mask(vector) for name, vector in axes.items()}
    d3 = []
    for composite, constituents in (
        ("nas_activity_sum", NAS_ACTIVITY),
        ("saf_activity_sum", SAF_ACTIVITY),
    ):
        any_constituent = np.zeros(universe.size, dtype=bool)
        for name in constituents:
            any_constituent |= masks[name]
        d3.append({
            "composite": composite,
            "constituents": list(constituents),
            "is_a_gate": False,
            "genes_bh_for_the_composite": int(masks[composite].sum()),
            "genes_bh_for_the_composite_and_no_constituent_alone": int(
                (masks[composite] & ~any_constituent).sum()
            ),
            "genes_bh_for_some_constituent_and_not_the_composite": int(
                (any_constituent & ~masks[composite]).sum()
            ),
            "genes_bh_for_any_constituent": int(any_constituent.sum()),
            "tie_structure_of_the_composite": ties[composite],
            "tie_structure_of_the_constituents": [ties[name] for name in constituents],
            "power_caveat": (
                "the composite has more distinct values and fewer tied pairs "
                "than any single constituent, so it has more power and a "
                "composite-only set can be a power difference rather than "
                "emergent signal. Never a gate."
            ),
        })
    print("d3 diagnostic done", flush=True)

    # PC1 sensitivity, never a check.
    def pc1_score(names: Sequence[str]) -> np.ndarray:
        block = np.column_stack([unit_ranks(axes[name]) for name in names])
        _, _, right = np.linalg.svd(block, full_matrices=False)
        loading = right[0]
        if loading.sum() < 0:
            loading = -loading
        return block @ loading

    pc1 = {}
    for composite, constituents in (
        ("nas_activity_sum", NAS_ACTIVITY),
        ("saf_activity_sum", SAF_ACTIVITY),
    ):
        score = pc1_score(constituents)
        pc1[composite] = {
            "is_a_gate": False,
            "spearman_with_the_crn_sum": spearman_correlation(
                score.tolist(), axes[composite].tolist()
            ),
            "marginal_genes_bh_below_0_05": marginal_family(
                unit_genes, score, axis=f"{composite}_pc1", n_participants=n
            )["genes_bh_below_0_05"],
            "crn_sum_marginal_genes_bh_below_0_05": marginal[composite][
                "genes_bh_below_0_05"
            ],
        }
    print("pc1 sensitivity done", flush=True)

    # Participant bootstrap of every criterion count, index-aligned so the two
    # directions of a criterion are paired on the same resample.
    bootstrap = bootstrap_counts(
        realised_values,
        axes,
        D1_FAMILIES + D2_FAMILIES,
        n_resamples=arguments.bootstrap,
        seed=arguments.seed + 3,
    )
    print("participant bootstrap done", flush=True)

    concentration = gain_concentration(
        [float(primary[f"{a}|{b}"]["genes_bh_below_0_05"]) for a, b in D1_FAMILIES],
        [
            float(primary[f"{a}|{b}"]["count_null"]["null_percentile_95_count"])
            for a, b in D1_FAMILIES
        ],
        ["nas_activity_sum|fibrosis", "fibrosis|nas_activity_sum"],
    )

    payload = {
        "schema_version": "masld-bench-axis-count-result-v1",
        "prespec_id": prespec["prespec_id"],
        "outcome": outcome,
        "grouping": grouping,
        "a_one_axis_answer_is_a_result_not_a_failure": True,
        "honest_provenance": prespec["honest_provenance"],
        "determination_audit_from_the_freeze": prespec[
            "determination_audit_run_before_the_freeze"
        ],
        "cohort": {
            "dataset_id": prespec["cohort"]["dataset_id"],
            "endpoint_rows": n,
            "distinct_participants": len(set(participants)),
            "unit_the_arithmetic_uses": "participant",
        },
        "tie_structure": ties,
        "criteria": {
            "d1_activity_and_fibrosis_are_independent": {
                "decisive": True,
                "threshold": FROZEN_THRESHOLDS["d1_activity_and_fibrosis_are_independent"],
                "gate": d1,
                "directions": [primary[f"{a}|{b}"] for a, b in D1_FAMILIES],
                "global_null_beside_the_primary": [
                    global_null[f"{a}|{b}"] for a, b in D1_FAMILIES
                ],
            },
            "d2_steatosis_and_saf_activity_are_independent": {
                "decisive": False,
                "threshold": FROZEN_THRESHOLDS[
                    "d2_steatosis_and_saf_activity_are_independent"
                ],
                "gate": d2,
                "directions": [primary[f"{a}|{b}"] for a, b in D2_FAMILIES],
                "global_null_beside_the_primary": [
                    global_null[f"{a}|{b}"] for a, b in D2_FAMILIES
                ],
                "recorded_limitation": prespec["criteria"][
                    "d2_steatosis_and_saf_activity_are_independent"
                ]["recorded_limitation"],
            },
        },
        "diagnostics_never_gates": {
            "d3_does_a_composite_earn_its_place": d3,
            "marginal_counts_for_every_axis": marginal,
            "pc1_composite_sensitivity": pc1,
            "participant_bootstrap_of_every_count": bootstrap,
            "label_level_context": label_context,
            "gain_concentration": {
                "applicable": concentration.applicable,
                "reason": concentration.reason,
                "n_strata": concentration.n_strata,
                "min_strata": concentration.min_strata,
                "supports_pass_fail_verdict": concentration.supports_pass_fail_verdict,
                "why_it_was_called_anyway": (
                    "so the not-applicable verdict is a measured return from "
                    "the package rather than an assertion in prose. There is "
                    "one study, one platform and one processing run here, so "
                    "there is no pooled gain across strata for it to decompose."
                ),
            },
            "paired_cluster_bootstrap_was_not_used": (
                "it reduces within cluster to one paired observation per "
                "biological unit, and a family-level gene count has no "
                "per-participant decomposition to reduce. The participant "
                "bootstrap above resamples participants directly and pairs the "
                "two directions on identical resample indices, which is the "
                "same guarantee for this statistic."
            ),
        },
        "method": prespec["method"],
        "power": {
            **prespec["power"],
            "measured_familywise_floors": {
                key: entry["familywise_floor"]["detectable_abs_partial_r_p95"]
                for key, entry in primary.items()
            },
        },
        "strata": {
            "recorded_sex": {
                "counts": dict(sorted(Counter(sexes).items())),
                "is_the_primary_reference": True,
                "why": (
                    "Stage 0 measured ballooning +0.353 at p=0.005 and lobular "
                    "inflammation +0.281 at p=0.010 against recorded_sex, so it "
                    "carries between-stratum structure here."
                ),
            },
            "outer_fold": {
                "counts": dict(sorted(Counter(folds.values()).items())),
                "is_a_genuine_nuisance_stratum": False,
                "why": (
                    "the folds are balanced on histology stage by construction, "
                    "so permuting within one would hold fixed part of the very "
                    "outcome variation the criteria ask about. Stage 0 ruled it "
                    "out and this run does not revisit that."
                ),
            },
        },
        "reuse_validation": {
            "rank_agreement_with_masld_bench": rank_agreement,
            "association_agreement_with_masld_bench": marginal_check,
            "bh_agreement_with_masld_bench": bh_agreement,
        },
        "controls": {
            "seed": arguments.seed,
            "count_null_draws": arguments.count_null_draws,
            "label_permutation_draws": arguments.label_draws,
            "bootstrap_resamples": arguments.bootstrap,
            "bootstrap_resampling_unit": "participant",
            "lobular_necrosis_never_entered_a_composite_or_a_criterion": True,
            "no_model_was_fitted": True,
            "nulls_measured_from_the_realised_residuals": True,
        },
        "claim_boundary": prespec["claim_boundary"],
    }

    arguments.output.mkdir(mode=0o750, parents=True)
    (arguments.output / "axis_count.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(arguments.output, {
        "artifact_class": "gse267145_axis_count_result",
        "prespec_id": prespec["prespec_id"],
        "outcome": outcome,
        "d1_passed": d1["passed"],
        "d2_passed": d2["passed"],
        "model_fitted": False,
        "status": "passed",
    })
    verify_frozen_tree(arguments.output)

    print(json.dumps({
        "outcome": outcome,
        "grouping": grouping,
        "d1": {"verdict": d1["verdict"],
               "counts": {entry["condition"]: entry["genes_bh_below_0_05"]
                          for entry in d1["conditions"]},
               "floors": {entry["condition"]: entry["count_null_percentile_95"]
                          for entry in d1["conditions"]}},
        "d2": {"verdict": d2["verdict"],
               "counts": {entry["condition"]: entry["genes_bh_below_0_05"]
                          for entry in d2["conditions"]},
               "floors": {entry["condition"]: entry["count_null_percentile_95"]
                          for entry in d2["conditions"]}},
        "realised_gene_universe": int(universe.size),
    }, indent=2))
    return 0


def bootstrap_counts(
    raw_genes: np.ndarray,
    axes: dict[str, np.ndarray],
    families: Sequence[tuple[str, str]],
    *,
    n_resamples: int,
    seed: int,
) -> dict[str, object]:
    """Participant bootstrap of every criterion count, on shared resamples.

    One resample of participants produces every family's count at once, so the
    intervals are comparable and a difference between two directions has a
    paired interval rather than two marginal ones. Ranks are recomputed inside
    each resample; genes that go constant there are dropped for that resample
    only and the surviving count is reported.
    """

    generator = np.random.default_rng(seed)
    n = raw_genes.shape[0]
    draws: dict[str, list[int]] = {f"{a}|{b}": [] for a, b in families}
    surviving: list[int] = []
    skipped = 0
    for _ in range(n_resamples):
        indices = generator.integers(0, n, size=n)
        block = raw_genes[indices, :]
        keep = block.min(axis=0) != block.max(axis=0)
        resampled = {name: vector[indices] for name, vector in axes.items()}
        if not keep.any() or any(
            vector.min() == vector.max() for vector in resampled.values()
        ):
            skipped += 1
            continue
        surviving.append(int(keep.sum()))
        ranks = average_ranks(block[:, keep])
        ranks -= ranks.mean(axis=0, keepdims=True)
        for exposure, covariate in families:
            unit_covariate = unit_ranks(resampled[covariate])
            residual_genes = residualize(ranks, unit_covariate)
            scaled, norms = unit_columns(residual_genes)
            usable = norms > 0.0
            residual_exposure = residualize(
                centred_ranks(resampled[exposure]), unit_covariate
            )
            norm = float(np.sqrt((residual_exposure**2).sum()))
            if norm == 0.0 or not usable.any():
                continue
            vector = scaled[:, usable].T @ (residual_exposure / norm)
            draws[f"{exposure}|{covariate}"].append(
                bh_count_from_abs_r(
                    np.abs(vector), bh_critical_abs_r(int(vector.size), n - 3)
                )
            )
    summary: dict[str, object] = {
        "resampling_unit": "participant",
        "n_resamples": n_resamples,
        "n_skipped_resamples": skipped,
        "median_genes_non_constant_in_the_resample": (
            float(np.median(surviving)) if surviving else float("nan")
        ),
        "seed": seed,
        "per_family": {},
        "paired_differences": {},
        "draws_are_index_aligned_across_families": True,
    }
    for key, counts in draws.items():
        array = np.asarray(counts, dtype=float)
        summary["per_family"][key] = {
            "n_usable_resamples": int(array.size),
            "median_count": float(np.median(array)) if array.size else None,
            "percentile_2_5": float(np.percentile(array, 2.5)) if array.size else None,
            "percentile_97_5": float(np.percentile(array, 97.5)) if array.size else None,
            "fraction_of_resamples_with_zero": (
                float(np.mean(array == 0.0)) if array.size else None
            ),
        }
    for left, right in (
        (f"{D1_FAMILIES[0][0]}|{D1_FAMILIES[0][1]}",
         f"{D1_FAMILIES[1][0]}|{D1_FAMILIES[1][1]}"),
        (f"{D2_FAMILIES[0][0]}|{D2_FAMILIES[0][1]}",
         f"{D2_FAMILIES[1][0]}|{D2_FAMILIES[1][1]}"),
    ):
        a = np.asarray(draws[left], dtype=float)
        b = np.asarray(draws[right], dtype=float)
        if a.size and a.size == b.size:
            difference = a - b
            summary["paired_differences"][f"{left} minus {right}"] = {
                "median": float(np.median(difference)),
                "percentile_2_5": float(np.percentile(difference, 2.5)),
                "percentile_97_5": float(np.percentile(difference, 97.5)),
                "why_paired": (
                    "the two directions are computed on identical resample "
                    "indices, so overlapping marginal intervals would not be "
                    "evidence of no difference"
                ),
            }
    return summary


if __name__ == "__main__":
    raise SystemExit(main())
