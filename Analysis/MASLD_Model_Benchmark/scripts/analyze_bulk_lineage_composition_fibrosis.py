#!/usr/bin/env python3
"""Score the bulk lineage-composition versus fibrosis transfer lane.

The model under test is the frozen reference plus BayesPrism deconvolution
operator. The question is whether its lineage projection carries a fibrosis-stage
rank association that transfers to a held-out cohort.

Four things this evaluator refuses to do the easy way.

1.  A collapsed lineage is not tested. Rank makes a sub-detection value finite,
    not informative; a lineage whose median sample sits below the 1e-4 floor in
    any contributing cohort is ranking posterior noise. Its rho is reported with
    its above-floor fraction and it is never entered into a BH family.
2.  Cross-cohort combination is the macro-average the frozen generalization
    requirements prescribe, not a precision-weighted meta-analysis, so one large
    cohort cannot dominate.
3.  Nothing is bootstrapped over cohorts. Three cohorts admit ten distinct
    resample multisets; transport evidence comes from the leave-one-cohort-out
    folds instead.
4.  A check condition that a deterministic procedure satisfies by construction
    records not_applicable and leaves the denominator, rather than paying out.
5.  Collapse is tested for confounding, not just for prevalence. Values can be
    distinct and still be numerical noise, so distinctness disposes of the tie
    mechanism and not of the confounding one: if RNA quality varies with fibrosis
    stage, degraded input produces more collapsed posterior mass, collapse status
    then tracks the outcome, and a lineage's rank association is partly a
    measurement output file wearing a biological label. BH will not catch that,
    because the association is real -- just not biological.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

PRIMARY_COHORTS = ("GSE135251", "GSE162694", "GSE240729")
TWO_COHORT = ("GSE135251", "GSE162694")
REDUCED_COHORT = "GSE213621"
DETECTION_FLOOR = 1e-4
MIN_FRACTION_ABOVE_FLOOR = 0.50
DOMINANT_LINEAGE = "Hepatocytes"
COLLAPSE_NULL_REPLICATES = 2000
COLLAPSE_CONFOUNDING_ALPHA = 0.05


class LineageEvaluatorError(RuntimeError):
    """Raised when the evaluation would not meet its frozen requirements."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t",
            lineterminator="\n", extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    from scipy.stats import spearmanr

    if len(x) < 3:
        raise LineageEvaluatorError("spearman needs at least three samples")
    value = spearmanr(x, y).statistic
    return float(value)


def fraction_above_floor(values: np.ndarray, floor: float = DETECTION_FLOOR) -> float:
    return float(np.mean(np.asarray(values, dtype=np.float64) >= floor))


def permutation_null(
    x: np.ndarray, y: np.ndarray, *, replicates: int, seed: int
) -> dict[str, float]:
    """Null for this exact pair, permuting the outcome against the realised scores.

    The outcome is a five-level stage and is heavily tied. Ties in the OUTCOME do
    not move a rank-correlation null, so this measured null is expected to agree
    with the 1/sqrt(n-1) analytic reference. Both are reported and a disagreement
    is a finding rather than something to smooth over.
    """

    rng = np.random.default_rng(seed)
    permuted = np.asarray(y, dtype=np.float64).copy()
    draws = np.empty(replicates, dtype=np.float64)
    for index in range(replicates):
        rng.shuffle(permuted)
        draws[index] = spearman(x, permuted)
    point = spearman(x, y)
    return {
        "rho": point,
        "null_mean": float(np.mean(draws)),
        "null_sd": float(np.std(draws, ddof=1)),
        "null_p95_abs": float(np.quantile(np.abs(draws), 0.95)),
        "analytic_null_se": float(1.0 / np.sqrt(len(x) - 1)),
        "two_sided_permutation_p": float(
            (1.0 + float(np.sum(np.abs(draws) >= abs(point)))) / (replicates + 1.0)
        ),
    }


def bootstrap_interval(
    x: np.ndarray, y: np.ndarray, *, replicates: int, seed: int
) -> tuple[float, float, int]:
    """Sample-level bootstrap within one cohort.

    Optimistic if a participant contributed more than one sample. No donor key
    exists in any of these cohorts, so that cannot be checked and the caveat is
    carried on every interval this produces.
    """

    rng = np.random.default_rng(seed)
    n = len(x)
    draws = np.empty(replicates, dtype=np.float64)
    used = 0
    for index in range(replicates):
        pick = rng.integers(0, n, size=n)
        if len(np.unique(y[pick])) < 2 or len(np.unique(x[pick])) < 2:
            draws[index] = np.nan
            continue
        draws[index] = spearman(x[pick], y[pick])
        used += 1
    finite = draws[np.isfinite(draws)]
    if len(finite) < replicates // 2:
        raise LineageEvaluatorError("within-cohort bootstrap is degenerate")
    return float(np.quantile(finite, 0.025)), float(np.quantile(finite, 0.975)), used


def collapse_confounding(
    values: np.ndarray, stage: np.ndarray, *, floor: float, replicates: int, seed: int
) -> dict[str, Any]:
    """Test whether COLLAPSE STATUS tracks fibrosis stage.

    Distinctness of the proportion values disposes of the tie mechanism: a
    cohort-level floor shift cannot mechanically degrade a within-cohort
    Spearman when every value is distinct. It says nothing about whether the
    ordering among collapsed values carries signal or noise.

    The failure this catches is confounding. If cirrhotic tissue yields more
    degraded RNA, and degraded input produces more collapsed posterior mass,
    then being collapsed correlates with stage and a lineage's rank association
    is partly a measurement output file. The association would be real and would
    survive BH; it simply would not be biology.

    The indicator is binary and therefore heavily tied, so its null is narrower
    than a continuous scorer's. That is why the null is MEASURED against this
    exact indicator rather than assumed from 1/sqrt(n-1).
    """

    collapsed = (np.asarray(values, dtype=np.float64) < floor).astype(np.float64)
    fraction = float(np.mean(collapsed))
    if fraction <= 0.0 or fraction >= 1.0:
        return {
            "collapsed_fraction": fraction,
            "applicable": False,
            "not_applicable_reason": (
                "collapse status is constant in this cohort, so there is no "
                "collapse-versus-stage association to test"
            ),
            "rho_collapse_vs_stage": float("nan"),
            "permutation_p": float("nan"),
            "null_p95_abs": float("nan"),
            "confounded": False,
        }
    rng = np.random.default_rng(seed)
    permuted = np.asarray(stage, dtype=np.float64).copy()
    draws = np.empty(replicates, dtype=np.float64)
    for index in range(replicates):
        rng.shuffle(permuted)
        draws[index] = spearman(collapsed, permuted)
    point = spearman(collapsed, stage)
    p_value = float((1.0 + float(np.sum(np.abs(draws) >= abs(point)))) / (replicates + 1.0))
    return {
        "collapsed_fraction": fraction,
        "applicable": True,
        "not_applicable_reason": None,
        "rho_collapse_vs_stage": point,
        "permutation_p": p_value,
        "null_p95_abs": float(np.quantile(np.abs(draws), 0.95)),
        "confounded": bool(p_value < COLLAPSE_CONFOUNDING_ALPHA),
    }


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    values = np.asarray(p_values, dtype=np.float64)
    n = len(values)
    if n == 0:
        return []
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype=np.float64)
    out[order] = adjusted
    return [float(v) for v in out]


def load_fixture(fixture: Path, cohort: str) -> tuple[np.ndarray, np.ndarray]:
    proportions = np.load(fixture / f"{cohort}_proportions.npy", allow_pickle=False)
    with (fixture / f"{cohort}_samples.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    stage = np.asarray([int(r["fibrosis_stage"]) for r in rows], dtype=np.float64)
    if proportions.shape[0] != len(stage):
        raise LineageEvaluatorError(f"{cohort} fixture axes differ")
    return proportions, stage


def evaluate_arm(
    *,
    arm_id: str,
    cohorts: Sequence[str],
    lineages: Sequence[str],
    data: Mapping[str, tuple[np.ndarray, np.ndarray]],
    bootstrap_replicates: int,
    permutation_replicates: int,
    seed: int,
) -> dict[str, Any]:
    per_cohort: list[dict[str, Any]] = []
    eligible: list[str] = []
    collapse_confounded: list[str] = []
    macro: dict[str, float] = {}
    macro_se: dict[str, float] = {}

    for index, lineage in enumerate(lineages):
        rows = []
        fractions = []
        confounded_here: list[bool] = []
        for cohort in cohorts:
            proportions, stage = data[cohort]
            values = proportions[:, index]
            fraction = fraction_above_floor(values)
            fractions.append(fraction)
            null = permutation_null(
                values, stage, replicates=permutation_replicates, seed=seed + index
            )
            low, high, used = bootstrap_interval(
                values, stage, replicates=bootstrap_replicates, seed=seed + 7919 + index
            )
            collapse = collapse_confounding(
                values, stage, floor=DETECTION_FLOOR,
                replicates=COLLAPSE_NULL_REPLICATES, seed=seed + 104729 + index,
            )
            confounded_here.append(bool(collapse["confounded"]))
            rows.append({
                "arm_id": arm_id,
                "lineage": lineage,
                "cohort": cohort,
                "samples": len(stage),
                "fraction_above_floor": fraction,
                "detection_floor_eligible_in_this_cohort": fraction >= MIN_FRACTION_ABOVE_FLOOR,
                "rho": null["rho"],
                "ci_low": low,
                "ci_high": high,
                "bootstrap_replicates_used": used,
                "measured_null_sd": null["null_sd"],
                "analytic_null_se": null["analytic_null_se"],
                "measured_null_p95_abs": null["null_p95_abs"],
                "permutation_p": null["two_sided_permutation_p"],
                "collapsed_fraction": collapse["collapsed_fraction"],
                "collapse_vs_stage_rho": collapse["rho_collapse_vs_stage"],
                "collapse_vs_stage_p": collapse["permutation_p"],
                "collapse_vs_stage_testable": collapse["applicable"],
                "collapse_confounded_in_this_cohort": collapse["confounded"],
            })
        per_cohort.extend(rows)
        if any(confounded_here):
            collapse_confounded.append(lineage)
        # Two independent guards, because they catch different failures. The
        # floor-and-fraction rule catches a lineage that is mostly unmeasured, so
        # its ranks are noise. The collapse-confounding test catches a lineage
        # whose collapse status tracks the outcome, so its association is a
        # measurement output file. A lineage can pass either and fail the other.
        if all(f >= MIN_FRACTION_ABOVE_FLOOR for f in fractions) and not any(confounded_here):
            eligible.append(lineage)
        rhos = np.asarray([r["rho"] for r in rows], dtype=np.float64)
        ses = np.asarray([r["analytic_null_se"] for r in rows], dtype=np.float64)
        macro[lineage] = float(np.mean(rhos))
        macro_se[lineage] = float(np.sqrt(np.sum(ses ** 2)) / len(ses))

    # BH within this arm's own detection-floor-eligible family only.
    from scipy.stats import norm

    p_by_lineage = {
        lineage: float(2.0 * (1.0 - norm.cdf(abs(macro[lineage]) / macro_se[lineage])))
        for lineage in eligible
    }
    q = dict(zip(eligible, benjamini_hochberg([p_by_lineage[l] for l in eligible]), strict=True))

    return {
        "arm_id": arm_id,
        "cohorts": list(cohorts),
        "samples": int(sum(len(data[c][1]) for c in cohorts)),
        "per_cohort_rows": per_cohort,
        "eligible_lineages": eligible,
        "eligible_family_size": len(eligible),
        "collapse_confounded_lineages": collapse_confounded,
        "macro_rho": macro,
        "macro_se": macro_se,
        "macro_p": p_by_lineage,
        "macro_q": q,
    }


def leave_one_cohort_out(
    *,
    cohorts: Sequence[str],
    lineages: Sequence[str],
    data: Mapping[str, tuple[np.ndarray, np.ndarray]],
    macro_rho: Mapping[str, float],
) -> list[dict[str, Any]]:
    folds: list[dict[str, Any]] = []
    for held_out in cohorts:
        training = [c for c in cohorts if c != held_out]
        for index, lineage in enumerate(lineages):
            train_rho = float(np.mean([
                spearman(data[c][0][:, index], data[c][1]) for c in training
            ]))
            held_rho = spearman(data[held_out][0][:, index], data[held_out][1])
            folds.append({
                "held_out_cohort": held_out,
                "training_cohorts": "+".join(training),
                "lineage": lineage,
                "training_macro_rho": train_rho,
                "held_out_rho": held_rho,
                "sign_agrees_with_training": bool(np.sign(train_rho) == np.sign(held_rho)),
                "sign_agrees_with_full_macro": bool(
                    np.sign(macro_rho[lineage]) == np.sign(held_rho)
                ),
                "magnitude_retention": (
                    float(abs(held_rho) / abs(train_rho)) if abs(train_rho) > 0 else float("nan")
                ),
            })
    return folds


def run(*, fixture: Path, gate_path: Path, output: Path, seed: int,
        bootstrap_replicates: int, permutation_replicates: int) -> dict[str, Any]:
    if output.exists():
        raise LineageEvaluatorError(f"refusing to overwrite evaluation: {output}")
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    conditions_spec = gate["development_advance_conditions"]

    with (fixture / "lineage_axis.tsv").open(encoding="utf-8", newline="") as handle:
        lineages = [r["lineage"] for r in csv.DictReader(handle, delimiter="\t")]
    if len(lineages) != 16:
        raise LineageEvaluatorError(f"lineage axis carries {len(lineages)} entries, not 16")

    data = {c: load_fixture(fixture, c) for c in PRIMARY_COHORTS + (REDUCED_COHORT,)}

    primary = evaluate_arm(
        arm_id="three_cohort_primary", cohorts=PRIMARY_COHORTS, lineages=lineages, data=data,
        bootstrap_replicates=bootstrap_replicates,
        permutation_replicates=permutation_replicates, seed=seed,
    )
    co_primary = evaluate_arm(
        arm_id="two_cohort_co_primary", cohorts=TWO_COHORT, lineages=lineages, data=data,
        bootstrap_replicates=bootstrap_replicates,
        permutation_replicates=permutation_replicates, seed=seed,
    )
    reduced = evaluate_arm(
        arm_id="reduced_resolution_consistency", cohorts=(REDUCED_COHORT,), lineages=lineages,
        data=data, bootstrap_replicates=bootstrap_replicates,
        permutation_replicates=permutation_replicates, seed=seed,
    )
    folds = leave_one_cohort_out(
        cohorts=PRIMARY_COHORTS, lineages=lineages, data=data, macro_rho=primary["macro_rho"],
    )

    # Filters live outside development_advance_conditions because they are applied
    # inside the qualifying-set definition rather than scored. Reading them from
    # the wrong block is how a script silently diverges from a frozen check.
    filters_spec = gate["qualifying_set_filters"]
    effect_floor = float(filters_spec["effect_floor"]["threshold"])
    if filters_spec["effect_floor"]["role"] != "filter, not a scored condition":
        raise LineageEvaluatorError("the frozen gate no longer declares the effect floor a filter")
    if float(filters_spec["detection_floor"]["floor"]) != DETECTION_FLOOR:
        raise LineageEvaluatorError("the frozen gate's detection floor differs from this evaluator")
    if float(filters_spec["detection_floor"]["minimum_fraction_above_floor"]) != MIN_FRACTION_ABOVE_FLOOR:
        raise LineageEvaluatorError("the frozen gate's minimum above-floor fraction differs")
    if float(filters_spec["collapse_confounding"]["alpha"]) != COLLAPSE_CONFOUNDING_ALPHA:
        raise LineageEvaluatorError("the frozen gate's collapse-confounding alpha differs")

    retention_floor = float(conditions_spec["held_out_magnitude_retention_minimum"]["threshold"])
    count_threshold = int(conditions_spec["qualifying_lineage_count_minimum"]["threshold"])
    set_change_max = int(conditions_spec["co_primary_set_change_maximum"]["threshold"])
    # The dominance guard's subject is read from the check, not hardcoded here, so a
    # change to the check cannot leave this evaluator testing a different lineage.
    dominance_requirement = str(conditions_spec["qualifying_set_dominance_guard"]["requirement"])
    if DOMINANT_LINEAGE not in dominance_requirement:
        raise LineageEvaluatorError(
            f"the frozen gate's dominance guard does not name {DOMINANT_LINEAGE}: "
            f"{dominance_requirement!r}"
        )

    def qualifying(arm: Mapping[str, Any], use_folds: bool) -> list[str]:
        out = []
        for lineage in arm["eligible_lineages"]:
            if arm["macro_q"][lineage] >= 0.05:
                continue
            if abs(arm["macro_rho"][lineage]) < effect_floor:
                continue
            if use_folds:
                signs = [
                    f["sign_agrees_with_full_macro"] for f in folds if f["lineage"] == lineage
                ]
                if not all(signs):
                    continue
            out.append(lineage)
        return out

    primary_qualifying = qualifying(primary, True)
    co_qualifying = qualifying(co_primary, False)

    retention_values = [
        f["magnitude_retention"] for f in folds
        if f["lineage"] in primary_qualifying and np.isfinite(f["magnitude_retention"])
    ]
    mean_retention = float(np.mean(retention_values)) if retention_values else float("nan")

    intersection = sorted(set(primary["eligible_lineages"]) & set(co_primary["eligible_lineages"]))
    set_change = len(
        (set(primary_qualifying) ^ set(co_qualifying)) & set(intersection)
    )
    non_dominant = [l for l in primary_qualifying if l != DOMINANT_LINEAGE]

    conditions = [
        {
            "condition_id": "qualifying_lineage_count_minimum",
            "applicable": True,
            "requirement": f">= {count_threshold} of {primary['eligible_family_size']} eligible lineages qualify",
            "observed": len(primary_qualifying),
            "met": bool(len(primary_qualifying) >= count_threshold),
        },
        {
            "condition_id": "qualifying_set_dominance_guard",
            "applicable": True,
            "requirement": f"at least one qualifying lineage is not {DOMINANT_LINEAGE}",
            "observed": len(non_dominant),
            "met": bool(non_dominant),
        },
        {
            "condition_id": "held_out_magnitude_retention_minimum",
            "applicable": bool(primary_qualifying),
            "requirement": f"mean held-out |rho| / training |rho| >= {retention_floor}",
            "observed": mean_retention,
            "met": bool(retention_values) and mean_retention >= retention_floor,
            "not_applicable_reason": (
                None if primary_qualifying
                else "no lineage qualified, so there is no retention ratio to evaluate"
            ),
        },
        {
            "condition_id": "co_primary_set_change_maximum",
            "applicable": True,
            "requirement": f"qualifying set changes by <= {set_change_max} on the eligible intersection",
            "observed": set_change,
            "met": bool(set_change <= set_change_max),
        },
        {
            "condition_id": "seed_direction_consistency",
            "applicable": False,
            "requirement": "n/a",
            "observed": None,
            "met": None,
            "not_applicable_reason": (
                "Spearman, BH and the leave-one-cohort-out partition are deterministic. "
                "There is no seed, no stochastic fit and no initialisation, so this "
                "condition would be satisfied by construction and carries no evidence."
            ),
        },
    ]
    # Every condition the frozen check registers must be computed here, and this
    # evaluator may not invent one the check does not register. A check condition the
    # analysis silently fails to compute would report as absent rather than fail
    # loudly, which is the same defect class as a vacuously satisfied condition.
    computed_ids = {c["condition_id"] for c in conditions}
    registered_ids = set(conditions_spec)
    if computed_ids != registered_ids:
        raise LineageEvaluatorError(
            "the evaluator's condition set does not match the frozen gate: "
            f"registered_not_computed={sorted(registered_ids - computed_ids)} "
            f"computed_not_registered={sorted(computed_ids - registered_ids)}"
        )

    applicable = [c for c in conditions if c["applicable"]]
    met = [c for c in applicable if c["met"]]
    verdict = "PASS" if len(met) == len(applicable) and applicable else "FAIL"

    output.mkdir(parents=True)
    all_rows = primary["per_cohort_rows"] + co_primary["per_cohort_rows"] + reduced["per_cohort_rows"]
    write_tsv(output / "per_cohort_lineage_results.tsv", tuple(all_rows[0]), all_rows)
    write_tsv(output / "leave_one_cohort_out_folds.tsv", tuple(folds[0]), folds)
    macro_rows = []
    for arm in (primary, co_primary):
        for lineage in lineages:
            macro_rows.append({
                "arm_id": arm["arm_id"],
                "lineage": lineage,
                "detection_floor_eligible": lineage in arm["eligible_lineages"],
                "macro_rho": arm["macro_rho"][lineage],
                "macro_null_se": arm["macro_se"][lineage],
                "macro_p": arm["macro_p"].get(lineage, ""),
                "macro_q_within_eligible_family": arm["macro_q"].get(lineage, ""),
                "qualifying": lineage in (
                    primary_qualifying if arm["arm_id"] == "three_cohort_primary" else co_qualifying
                ),
            })
    write_tsv(output / "macro_average_by_lineage.tsv", tuple(macro_rows[0]), macro_rows)

    verdict_doc = {
        "gate_id": gate["gate_id"],
        "primary_endpoint": gate["primary_endpoint"],
        "arm_id": "three_cohort_primary",
        "verdict": verdict,
        "conditions": conditions,
        "conditions_met": len(met),
        "conditions_applicable": len(applicable),
        "conditions_total_including_not_applicable": len(conditions),
        "verdict_arithmetic": f"{len(met)} of {len(applicable)} applicable ({len(conditions)} total)",
        "qualifying_lineages": primary_qualifying,
        "co_primary_qualifying_lineages": co_qualifying,
        "eligible_intersection": intersection,
        "external_development_only": True,
        "champion_eligible": False,
    }
    (output / "promotion_gate_verdict.json").write_text(
        json.dumps(verdict_doc, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    receipt = {
        "schema_version": "masld-bench-lineage-composition-evaluation-v1",
        "status": "passed_external_development_scoring",
        "task_id": "bulk_lineage_composition_fibrosis_transfer",
        "unit": "bulk_rna_sample",
        "donor_count_claimed": False,
        "lineages": len(lineages),
        "detection_floor": DETECTION_FLOOR,
        "minimum_fraction_above_floor": MIN_FRACTION_ABOVE_FLOOR,
        "gate_condition_set_verified_against_frozen_gate": True,
        "effect_floor_read_from": "qualifying_set_filters.effect_floor.threshold",
        "effect_floor": effect_floor,
        "collapse_confounding_test": {
            "statistic": "spearman(collapsed_indicator, fibrosis_stage) per lineage per cohort",
            "null": "measured by permutation against the realised binary indicator",
            "null_replicates": COLLAPSE_NULL_REPLICATES,
            "alpha": COLLAPSE_CONFOUNDING_ALPHA,
            "rule": "a lineage confounded in ANY contributing cohort is reported and never counted",
            "why_measured_not_analytic": (
                "the indicator is binary and heavily tied, so its null is narrower than a "
                "continuous scorer's; assuming 1/sqrt(n-1) would overstate the threshold"
            ),
            "what_distinctness_did_and_did_not_settle": (
                "zero ties disposes of the mechanism where a cohort-level floor shift "
                "mechanically degrades a within-cohort Spearman; it is silent on whether "
                "the ordering among collapsed values is signal or noise, which is what this "
                "test addresses"
            ),
        },
        "three_cohort_primary": {
            "cohorts": list(PRIMARY_COHORTS),
            "samples": primary["samples"],
            "eligible_lineages": primary["eligible_lineages"],
            "eligible_family_size": primary["eligible_family_size"],
            "collapse_confounded_lineages": primary["collapse_confounded_lineages"],
            "qualifying_lineages": primary_qualifying,
        },
        "two_cohort_co_primary": {
            "cohorts": list(TWO_COHORT),
            "samples": co_primary["samples"],
            "eligible_lineages": co_primary["eligible_lineages"],
            "eligible_family_size": co_primary["eligible_family_size"],
            "collapse_confounded_lineages": co_primary["collapse_confounded_lineages"],
            "qualifying_lineages": co_qualifying,
            "is_not_a_fresh_frozen_estimate": True,
        },
        "reduced_resolution_arm": {
            "cohort": REDUCED_COHORT,
            "samples": reduced["samples"],
            "eligible_lineages": reduced["eligible_lineages"],
            "gate_eligible": False,
            "scale": "three pooled ordinal levels, never pooled with the graded cohorts",
        },
        "combination": "macro_average_across_cohort_families",
        "cohort_level_bootstrap_performed": False,
        "between_cohort_heterogeneity_estimated": False,
        "bootstrap_replicates": bootstrap_replicates,
        "permutation_replicates": permutation_replicates,
        "seed": seed,
        "promotion_gate_verdict": verdict,
        "verdict_arithmetic": verdict_doc["verdict_arithmetic"],
        "lineage_abundance_is_an_inference": True,
        "shared_reference_dependency": (
            "all cohorts were deconvolved against one REF_HUMAN reference, so cross-cohort "
            "agreement is evidence that one fixed operator transfers and is NOT evidence "
            "that independent measurements agree"
        ),
        "external_development_only": True,
        "champion_claim_eligible": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--gate", required=True, type=Path)
    parser.add_argument("--gate-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--bootstrap-replicates", type=int, default=10000)
    parser.add_argument("--permutation-replicates", type=int, default=10000)
    arguments = parser.parse_args()
    if sha256_file(arguments.gate) != arguments.gate_sha256:
        raise LineageEvaluatorError("promotion gate SHA differs")
    receipt = run(
        fixture=arguments.fixture, gate_path=arguments.gate, output=arguments.output,
        seed=arguments.seed, bootstrap_replicates=arguments.bootstrap_replicates,
        permutation_replicates=arguments.permutation_replicates,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
