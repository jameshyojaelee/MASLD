"""Answer the frozen Arm B axis-count question, prespecification proved first.

Stage 0c of the MASLD showcase model, job 2 of 2. The wrapper verifies the
frozen tree with ``sha256sum --check --strict`` before this runs, and this
script re-reads the criteria out of that frozen file rather than restating them.

Every statistic is imported from ``evaluate_gse267145_axis_count.py`` rather
than reimplemented, so Arm A and Arm B are compared on the same instrument. The
only deliberate differences from Stage 0b are recorded in the frozen document:
the primary permutation is global because GSE135251 deposits no biological
stratum, and the participant bootstrap is replaced by a leave-one-out jackknife
because resampling with replacement was shown to be invalid for a BH count.

E1 decides whether both directions survive. E2 classifies any empty direction
against its own measured familywise floor -- underpowered below it, genuinely
not independent at or above it -- which is the distinction that settled Stage
0b's interpretation and is stated in advance here rather than found again.
"""

from __future__ import annotations

import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class ArmBError(RuntimeError):
    """Raised when the substrate cannot support the frozen question."""


ACTIVITY = "nas_score"
FIBROSIS = "fibrosis_stage"
FAMILIES = ((ACTIVITY, FIBROSIS), (FIBROSIS, ACTIVITY))

FROZEN_THRESHOLDS = {
    "e1_activity_and_fibrosis_are_independent": (
        "in each direction, BH 0.05 count > 0 and above the 95th percentile of "
        "the matched permutation count null"
    ),
    "e2_an_empty_direction_is_classified_against_the_measured_floor": (
        "an empty direction is underpowered if its observed max |partial r| is "
        "below the measured familywise floor p95, and not_independent if it "
        "reaches or clears that floor"
    ),
}


def _import(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ArmBError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_prespecification(directory: Path) -> dict:
    verify_frozen_tree(directory)
    payload = json.loads(
        (directory / "arm_b_axis_count_prespecification.json").read_text(
            encoding="utf-8"
        )
    )
    criteria = payload["criteria"]
    if sorted(criteria) != sorted(FROZEN_THRESHOLDS):
        raise ArmBError("the frozen criteria are not the two expected ones")
    for name, expected in FROZEN_THRESHOLDS.items():
        if criteria[name]["threshold"] != expected:
            raise ArmBError(f"{name} threshold drifted since the freeze")
    if criteria["e1_activity_and_fibrosis_are_independent"].get("decisive") is not True:
        raise ArmBError("e1 lost its decisive mark since the freeze")
    if payload["honest_provenance"][
        "criteria_fixed_before_any_arm_b_number_was_computed"
    ] is not True:
        raise ArmBError("the frozen provenance does not assert an untouched Arm B")
    return payload


def classify_empty_direction(result: dict) -> dict[str, object]:
    """E2. An empty count means one of two very different things.

    Below the measured familywise floor the best gene in the family never had a
    chance of being detected, so the emptiness is a statement about the
    substrate. At or above the floor, genes at that level were detectable and
    none survived adjustment, so the emptiness is evidence.
    """

    observed_max = float(result["observed_max_abs_partial_r"])
    floor = float(result["familywise_floor"]["detectable_abs_partial_r_p95"])
    applicable = int(result["genes_bh_below_0_05"]) == 0
    if not applicable:
        return {
            "direction": f"{result['exposure']} | {result['adjusted_for']}",
            "applicable": False,
            "why_not_applicable": "the direction is not empty, so there is nothing to disambiguate",
            "classification": None,
        }
    underpowered = observed_max < floor
    return {
        "direction": f"{result['exposure']} | {result['adjusted_for']}",
        "applicable": True,
        "observed_max_abs_partial_r": observed_max,
        "measured_familywise_floor_p95": floor,
        "classification": "underpowered" if underpowered else "not_independent",
        "reading": (
            "the best gene in this family sits below the level at which any "
            "gene could have been detected, so the empty count is a statement "
            "about the substrate and the axis is indeterminate, never "
            "tested_negative"
            if underpowered
            else "genes at the detectable level existed and none survived "
            "adjustment, so the empty count is evidence rather than an absence "
            "of evidence"
        ),
    }


def jackknife_counts(
    analysis,
    raw_genes: np.ndarray,
    axes: dict[str, np.ndarray],
    families: Sequence[tuple[str, str]],
) -> dict[str, object]:
    """Leave-one-participant-out influence on each count.

    Stage 0b's participant bootstrap was invalid for this statistic: every
    observed count sat below its own bootstrap 2.5th percentile, because
    sampling with replacement duplicates participants whose identical
    expression and identical labels manufacture exact concordance blocks that
    inflate rank correlation. A jackknife introduces no duplicates.
    """

    n = raw_genes.shape[0]
    per_family: dict[str, list[int]] = {f"{a}|{b}": [] for a, b in families}
    for dropped in range(n):
        keep_rows = np.arange(n) != dropped
        block = raw_genes[keep_rows, :]
        keep = block.min(axis=0) != block.max(axis=0)
        ranks = analysis.average_ranks(block[:, keep])
        ranks -= ranks.mean(axis=0, keepdims=True)
        for exposure, covariate in families:
            values = {name: vector[keep_rows] for name, vector in axes.items()}
            if values[covariate].min() == values[covariate].max():
                continue
            unit_covariate = analysis.unit_ranks(values[covariate])
            scaled, norms = analysis.unit_columns(
                analysis.residualize(ranks, unit_covariate)
            )
            usable = norms > 0.0
            residual = analysis.residualize(
                analysis.centred_ranks(values[exposure]), unit_covariate
            )
            norm = float(np.sqrt((residual**2).sum()))
            if norm == 0.0 or not usable.any():
                continue
            vector = scaled[:, usable].T @ (residual / norm)
            per_family[f"{exposure}|{covariate}"].append(
                analysis.bh_count_from_abs_r(
                    np.abs(vector),
                    analysis.bh_critical_abs_r(int(vector.size), n - 1 - 3),
                )
            )
    summary: dict[str, object] = {
        "method": "leave_one_participant_out",
        "why_not_a_bootstrap": (
            "resampling with replacement duplicates participants, and a "
            "duplicated participant carries identical expression and identical "
            "labels, manufacturing exact concordance blocks that inflate rank "
            "correlation. In Stage 0b every observed count fell below its own "
            "bootstrap 2.5th percentile, which is the signature of that bias."
        ),
        "n_leave_one_out_fits": n,
        "per_family": {},
    }
    for key, counts in per_family.items():
        array = np.asarray(counts, dtype=float)
        summary["per_family"][key] = {
            "n_fits": int(array.size),
            "median": float(np.median(array)) if array.size else None,
            "minimum": float(array.min()) if array.size else None,
            "maximum": float(array.max()) if array.size else None,
            "fraction_of_fits_with_zero": (
                float(np.mean(array == 0.0)) if array.size else None
            ),
        }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prespecification", type=Path, required=True)
    parser.add_argument("--analysis-script", type=Path, required=True)
    parser.add_argument("--addendum-script", type=Path, required=True)
    parser.add_argument("--endpoints", type=Path, required=True)
    parser.add_argument("--molecular", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count-null-draws", type=int, default=20000)
    parser.add_argument("--null-block", type=int, default=250)
    parser.add_argument("--splits", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260827)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise ArmBError("refusing to overwrite an Arm B axis-count result")

    analysis = _import(arguments.analysis_script, "arm_a_analysis")
    addendum = _import(arguments.addendum_script, "arm_a_addendum")
    prespec = load_prespecification(arguments.prespecification)
    print("verified the frozen prespecification", flush=True)

    rows = analysis.read_table(arguments.endpoints)
    participants = [row["participant_id"] for row in rows]
    n = len(rows)
    if n != prespec["cohort"]["n_participants"]:
        raise ArmBError("the participant count differs from the freeze")
    if len(set(participants)) != n:
        raise ArmBError("participant_id is not unique")

    axes = {
        ACTIVITY: np.asarray([int(row[ACTIVITY]) for row in rows], dtype=float),
        FIBROSIS: np.asarray([int(row[FIBROSIS]) for row in rows], dtype=float),
    }
    if int(axes[FIBROSIS].max()) != 4:
        raise ArmBError("Arm B fibrosis must realise stage 4")
    recorded = prespec["what_was_inspected_before_freezing"][
        "marginals_and_tie_structure_only"
    ]
    for name, vector in axes.items():
        observed = {str(k): v for k, v in sorted(Counter(vector.astype(int).tolist()).items())}
        if observed != recorded[name]["value_counts"]:
            raise ArmBError(f"{name} marginal differs from the frozen prespecification")

    axis_rows = analysis.read_table(arguments.molecular / "participant_axis.tsv")
    if [row["participant_id"] for row in axis_rows] != participants:
        raise ArmBError("the molecular participant axis is not in endpoint order")
    features = analysis.read_table(arguments.molecular / "rna_feature_axis.tsv")
    gene_ids = [row["stable_gene_id"] for row in features]
    if len(set(gene_ids)) != len(gene_ids):
        raise ArmBError("stable_gene_id is not unique on the feature axis")
    values = np.load(arguments.molecular / "rna_values.npy")
    if values.shape != (n, len(gene_ids)):
        raise ArmBError("the expression matrix does not match its axes")
    print(f"loaded expression {values.shape}", flush=True)

    ordered = np.sort(values, axis=0)
    realised = ((ordered[1:] != ordered[:-1]).sum(axis=0) + 1) >= 2
    del ordered
    universe = np.flatnonzero(realised)
    if universe.size == 0:
        raise ArmBError("no gene varies across participants")
    realised_values = values[:, universe]
    print(
        f"realised gene universe: {universe.size} of {len(gene_ids)} "
        f"({len(gene_ids) - universe.size} constant across all {n} participants)",
        flush=True,
    )

    generator = np.random.default_rng(arguments.seed)
    sample = sorted(
        generator.choice(universe.size, size=min(300, universe.size), replace=False).tolist()
    )
    rank_agreement = analysis.validate_rank_agreement(realised_values, sample)
    if not rank_agreement["agrees_with_masld_bench_average_ranks"]:
        raise ArmBError("the vectorized ranks disagree with the package")

    gene_ranks = analysis.average_ranks(realised_values)
    gene_ranks -= gene_ranks.mean(axis=0, keepdims=True)
    unit_genes, norms = analysis.unit_columns(gene_ranks)
    if np.any(norms == 0.0):
        raise ArmBError("a constant column survived the universe filter")

    marginal_check = []
    from masld_bench.evaluators.metrics import spearman_correlation
    for name in axes:
        vector = unit_genes.T @ analysis.unit_ranks(axes[name])
        worst = max(
            abs(
                spearman_correlation(
                    realised_values[:, column].tolist(), axes[name].tolist()
                )
                - float(vector[column])
            )
            for column in sample[:100]
        )
        marginal_check.append({
            "axis": name,
            "genes_checked": 100,
            "max_absolute_difference_vs_masld_bench_spearman": worst,
        })
    if max(e["max_absolute_difference_vs_masld_bench_spearman"] for e in marginal_check) > 1e-10:
        raise ArmBError("the vectorized association disagrees with the package")
    print("rank and association agreement validated", flush=True)

    directions = {}
    for offset, (exposure, covariate) in enumerate(FAMILIES):
        unit_covariate = analysis.unit_ranks(axes[covariate])
        scaled, residual_norms = analysis.unit_columns(
            analysis.residualize(gene_ranks, unit_covariate)
        )
        keep = residual_norms > 0.0
        result = analysis.partial_family(
            scaled[:, keep],
            analysis.residualize(analysis.centred_ranks(axes[exposure]), unit_covariate),
            unit_covariate,
            exposure=exposure,
            covariate=covariate,
            n_participants=n,
            n_draws=arguments.count_null_draws,
            block=arguments.null_block,
            seed=arguments.seed + 10 + offset,
            strata=None,
            null_name="global_residual_permutation",
        )
        result["realised_universe"] = {
            "features_on_axis": len(gene_ids),
            "non_constant_across_all_participants": int(universe.size),
            "dropped_as_collinear_with_the_covariate": int((~keep).sum()),
            "entered_this_family": int(keep.sum()),
            "re_derived_not_copied": True,
        }
        directions[f"{exposure}|{covariate}"] = result
        print(f"family done: {exposure}|{covariate}", flush=True)

    e1 = analysis.evaluate_gate(
        "e1_activity_and_fibrosis_are_independent",
        [analysis.direction_condition(directions[f"{a}|{b}"]) for a, b in FAMILIES],
    )
    e2 = analysis.evaluate_gate(
        "e2_an_empty_direction_is_classified_against_the_measured_floor",
        [
            {**classify_empty_direction(directions[f"{a}|{b}"]), "met": True}
            for a, b in FAMILIES
        ],
    )
    e2_by_direction = {
        entry["direction"]: entry for entry in e2["conditions"] if entry["applicable"]
    }

    survivors = [
        entry["condition"] for entry in e1["conditions"]
        if entry["applicable"] and entry["met"]
    ]
    if e1["verdict"] == "NO_APPLICABLE_CONDITIONS":
        outcome, meaning = "INDETERMINATE", "no applicable condition"
    elif e1["passed"]:
        outcome = "TWO_AXIS_ACTIVITY_FIBROSIS"
        meaning = "activity and fibrosis each carry genes the other does not explain"
    elif not survivors:
        outcome = "NO_AXIS_RESOLVED"
        meaning = "neither axis carries genes the other does not explain"
    else:
        survivor = survivors[0]
        other = f"{FIBROSIS} | {ACTIVITY}" if survivor.startswith(ACTIVITY) else f"{ACTIVITY} | {FIBROSIS}"
        empty = e2_by_direction.get(other)
        if empty is None:
            outcome, meaning = "INDETERMINATE", "the non-surviving direction is not empty and not met"
        elif survivor.startswith(ACTIVITY):
            outcome = (
                "ONE_AXIS_ACTIVITY_FIBROSIS_UNDERPOWERED"
                if empty["classification"] == "underpowered"
                else "ONE_AXIS_ACTIVITY_FIBROSIS_NOT_INDEPENDENT"
            )
            meaning = empty["reading"]
        else:
            outcome = (
                "ONE_AXIS_FIBROSIS_ACTIVITY_UNDERPOWERED"
                if empty["classification"] == "underpowered"
                else "ONE_AXIS_FIBROSIS_ACTIVITY_NOT_INDEPENDENT"
            )
            meaning = empty["reading"]
    print(f"outcome: {outcome}", flush=True)

    marginal = {
        name: analysis.marginal_family(unit_genes, vector, axis=name, n_participants=n)
        for name, vector in axes.items()
    }
    reliability = addendum.split_half_reliability(
        analysis, realised_values, axes, n_splits=arguments.splits, seed=arguments.seed + 1
    )
    arm_a = prespec["arm_a_prior_observations"][
        "split_half_gene_ordering_reliability_full_length"
    ]
    fibrosis_reliability = reliability[FIBROSIS]["spearman_brown_full_length_reliability"]
    print("reliability and marginals done", flush=True)

    jackknife = jackknife_counts(analysis, realised_values, axes, FAMILIES)
    print("jackknife done", flush=True)

    payload = {
        "schema_version": "masld-bench-arm-b-axis-count-result-v1",
        "prespec_id": prespec["prespec_id"],
        "outcome": outcome,
        "meaning": meaning,
        "a_one_axis_answer_is_a_result_not_a_failure": True,
        "honest_provenance": prespec["honest_provenance"],
        "arm_a_prior_observations": prespec["arm_a_prior_observations"],
        "fibrosis_scale_hazard": prespec["fibrosis_scale_hazard"],
        "cohort": {
            "dataset_id": prespec["cohort"]["dataset_id"],
            "endpoint_rows": n,
            "distinct_participants": len(set(participants)),
            "unit_the_arithmetic_uses": "participant",
        },
        "tie_structure": {
            name: analysis.tie_profile(vector, name) for name, vector in axes.items()
        },
        "criteria": {
            "e1_activity_and_fibrosis_are_independent": {
                "decisive": True,
                "threshold": FROZEN_THRESHOLDS["e1_activity_and_fibrosis_are_independent"],
                "gate": e1,
                "directions": [directions[f"{a}|{b}"] for a, b in FAMILIES],
            },
            "e2_an_empty_direction_is_classified_against_the_measured_floor": {
                "decisive": False,
                "is_a_disambiguation_rule_not_a_pass_fail_gate": True,
                "threshold": FROZEN_THRESHOLDS[
                    "e2_an_empty_direction_is_classified_against_the_measured_floor"
                ],
                "gate": e2,
                "note": (
                    "E2 never overrides E1. It qualifies the outcome name for a "
                    "direction E1 has already found empty."
                ),
            },
        },
        "diagnostics_never_gates": {
            "split_half_reliability_per_axis": {
                "is_a_gate": False,
                "per_axis": reliability,
                "arm_a_reference": arm_a,
                "prediction_recorded_before_the_measurement": prespec[
                    "diagnostics_never_gates"]["split_half_reliability_per_axis"][
                    "prediction_recorded_before_the_measurement"],
                "prediction_held": bool(fibrosis_reliability > arm_a["fibrosis"]),
                "arm_b_fibrosis_full_length_reliability": fibrosis_reliability,
                "arm_a_fibrosis_full_length_reliability": arm_a["fibrosis"],
            },
            "marginal_counts_for_both_axes": marginal,
            "participant_jackknife": jackknife,
        },
        "method": prespec["method"],
        "power": {
            **prespec["power"],
            "measured_familywise_floors": {
                key: entry["familywise_floor"]["detectable_abs_partial_r_p95"]
                for key, entry in directions.items()
            },
        },
        "reuse_validation": {
            "rank_agreement_with_masld_bench": rank_agreement,
            "association_agreement_with_masld_bench": marginal_check,
            "statistics_imported_from_stage_0b_not_reimplemented": True,
        },
        "controls": {
            **prespec["controls"],
            "seed": arguments.seed,
            "count_null_draws": arguments.count_null_draws,
            "split_half_repeats": arguments.splits,
            "primary_null_is_global_because_no_biological_stratum_is_deposited": True,
        },
        "claim_boundary": prespec["claim_boundary"],
    }

    arguments.output.mkdir(mode=0o750, parents=True)
    (arguments.output / "arm_b_axis_count.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(arguments.output, {
        "artifact_class": "gse135251_axis_count_result",
        "prespec_id": prespec["prespec_id"],
        "outcome": outcome,
        "e1_passed": e1["passed"],
        "model_fitted": False,
        "status": "passed",
    })
    verify_frozen_tree(arguments.output)

    print(json.dumps({
        "outcome": outcome,
        "e1_verdict": e1["verdict"],
        "counts": {
            entry["condition"]: entry["genes_bh_below_0_05"] for entry in e1["conditions"]
        },
        "floors": {
            key: round(entry["familywise_floor"]["detectable_abs_partial_r_p95"], 4)
            for key, entry in directions.items()
        },
        "observed_max": {
            key: round(entry["observed_max_abs_partial_r"], 4)
            for key, entry in directions.items()
        },
        "e2": {k: v["classification"] for k, v in e2_by_direction.items()},
        "arm_b_fibrosis_reliability": round(fibrosis_reliability, 4),
        "realised_gene_universe": int(universe.size),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
