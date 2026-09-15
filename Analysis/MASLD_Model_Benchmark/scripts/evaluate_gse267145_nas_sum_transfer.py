#!/usr/bin/env python3
"""Independently score the frozen GSE267145 activity-sum transfer.

Four commitments carried over from the frozen TaskSpec and check:

1.  Spearman is primary because it is invariant to the declared 1.688-point
    severity shift between the two cohorts.  Mean absolute error and the
    calibration slope and intercept are reported and are barred from advancing
    or blocking the lane, because they are not invariant to that shift and
    screening on them would fail every candidate for a property of the cohorts.
2.  The null is the permutation null of each model's own realised score vector,
    reported beside the analytic 1/sqrt(n-1) reference.  Outcome ties do not
    narrow a rank-correlation null and the 24 tied zeros are not used to argue
    for a tighter one.
3.  Intervals come from a participant bootstrap and the candidate-versus-
    baseline comparison is paired on the same resampled participants.
4.  A deterministic fit cannot earn the seed-direction condition.  When the
    minimum pairwise correlation across seed vectors exceeds the check's
    threshold the condition is recorded not_applicable, and the verdict reports
    met out of applicable so that three of three can never read as four of four.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

CLASSES_FREE = True
EXPECTED_PARTICIPANTS = 99
MANDATORY_BASELINES = (
    "gene_median_pca_ridge", "gene_median_ridge", "per_array_rank_pca_ridge",
    "per_array_rank_ridge", "training_mean_nas_sum",
)
PRIMARY_METRIC = "participant_spearman_nash_crn_component_sum"
METRIC_IDS = (
    PRIMARY_METRIC,
    "participant_mae_nash_crn_component_sum",
    "participant_calibration_slope_nash_crn_component_sum",
    "participant_calibration_intercept_nash_crn_component_sum",
)
NON_GATED_METRICS = METRIC_IDS[1:]


class NasEvaluatorError(RuntimeError):
    """Raised when predictions or outcomes does not meet the frozen evaluator requirements."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise NasEvaluatorError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _ranks(values: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    return rankdata(np.asarray(values, dtype=np.float64), method="average")


def spearman(observed: np.ndarray, predicted: np.ndarray) -> float:
    x, y = _ranks(observed), _ranks(predicted)
    x = x - x.mean()
    y = y - y.mean()
    nx, ny = float(np.linalg.norm(x)), float(np.linalg.norm(y))
    if nx <= 0 or ny <= 0:
        return 0.0
    return float(np.dot(x, y) / (nx * ny))


def mean_absolute_error(observed: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.mean(np.abs(np.asarray(predicted, dtype=np.float64) - observed)))


def calibration(observed: np.ndarray, predicted: np.ndarray) -> tuple[float, float]:
    """OLS of observed on predicted: slope 1 and intercept 0 is perfect."""

    x = np.asarray(predicted, dtype=np.float64)
    y = np.asarray(observed, dtype=np.float64)
    # np.var of an exactly constant vector returns floating-point residue near
    # 1e-31 rather than zero, which would divide a zero covariance by it and
    # report a calibration line for a predictor that has none.  Count distinct
    # values instead: it is exact and scale-free.
    if len(np.unique(x)) < 2:
        return float("nan"), float("nan")
    variance = float(np.var(x))
    if variance <= 0:
        return float("nan"), float("nan")
    slope = float(np.cov(x, y, ddof=0)[0, 1] / variance)
    return slope, float(np.mean(y) - slope * np.mean(x))


def metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    slope, intercept = calibration(observed, predicted)
    return {
        PRIMARY_METRIC: spearman(observed, predicted),
        "participant_mae_nash_crn_component_sum": mean_absolute_error(observed, predicted),
        "participant_calibration_slope_nash_crn_component_sum": slope,
        "participant_calibration_intercept_nash_crn_component_sum": intercept,
    }


def permutation_null(observed: np.ndarray, predicted: np.ndarray, *,
                     replicates: int, seed: int) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    score = _ranks(predicted)
    score = score - score.mean()
    score_norm = float(np.linalg.norm(score))
    labels = _ranks(observed)
    labels = labels - labels.mean()
    label_norm = float(np.linalg.norm(labels))
    n = len(observed)
    analytic = 1.0 / np.sqrt(n - 1)
    if score_norm <= 0 or label_norm <= 0:
        return {
            "null_mean": 0.0, "null_p95": 0.0, "null_sd": 0.0,
            "permutation_p_value": 1.0, "distinct_prediction_values": 1,
            "analytic_null_sd": float(analytic),
            "permutation_replicates": int(replicates),
            "p_value_resolution_floor": float(1.0 / (replicates + 1.0)),
        }
    working = labels.copy()
    draws = np.empty(replicates, dtype=np.float64)
    for index in range(replicates):
        rng.shuffle(working)
        draws[index] = float(np.dot(working, score) / (score_norm * label_norm))
    point = float(np.dot(labels, score) / (score_norm * label_norm))
    return {
        "null_mean": float(np.mean(draws)),
        "null_p95": float(np.quantile(draws, 0.95)),
        "null_sd": float(np.std(draws, ddof=1)),
        "analytic_null_sd": float(analytic),
        "permutation_p_value": float((1.0 + float(np.sum(draws >= point))) / (replicates + 1.0)),
        "permutation_replicates": int(replicates),
        "p_value_resolution_floor": float(1.0 / (replicates + 1.0)),
        "distinct_prediction_values": int(len(np.unique(np.round(predicted, 12)))),
    }


def load_predictions(path: Path, row_ids: Sequence[str]) -> tuple[str, np.ndarray]:
    fields, rows = read_tsv(path)
    if not {"row_id", "model_id", "predicted_nas_sum"} <= set(fields) or len(rows) != len(row_ids):
        raise NasEvaluatorError("prediction schema or row count differs")
    model_ids = {row["model_id"] for row in rows}
    if len(model_ids) != 1:
        raise NasEvaluatorError("prediction bundle contains multiple models")
    by_id = {row["row_id"]: row for row in rows}
    if len(by_id) != len(rows) or set(by_id) != set(row_ids):
        raise NasEvaluatorError("prediction row IDs differ")
    values = np.asarray([float(by_id[v]["predicted_nas_sum"]) for v in row_ids])
    if not np.all(np.isfinite(values)):
        raise NasEvaluatorError("predictions are not finite")
    return next(iter(model_ids)), values


def load_seed_scores(path: Path, row_ids: Sequence[str]) -> tuple[list[int], np.ndarray]:
    fields, rows = read_tsv(path)
    if not {"row_id", "model_seed", "predicted_nas_sum"} <= set(fields):
        raise NasEvaluatorError("seed-score schema differs")
    seeds = sorted({int(row["model_seed"]) for row in rows})
    table = {(row["row_id"], int(row["model_seed"])): float(row["predicted_nas_sum"])
             for row in rows}
    if len(table) != len(rows) or len(rows) != len(seeds) * len(row_ids):
        raise NasEvaluatorError("seed-score table is not rectangular")
    return seeds, np.asarray([[table[(r, s)] for r in row_ids] for s in seeds])


def min_pairwise_correlation(stack: np.ndarray) -> float:
    if stack.shape[0] < 2:
        return 1.0
    centred = stack - stack.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centred, axis=1, keepdims=True)
    if np.any(norms <= 0):
        return 1.0
    unit = centred / norms
    corr = unit @ unit.T
    return float(np.min(corr[np.triu_indices(stack.shape[0], k=1)]))


def load_outcomes(path: Path) -> tuple[list[str], np.ndarray]:
    fields, rows = read_tsv(path)
    if not {"row_id", "nash_crn_component_sum"} <= set(fields):
        raise NasEvaluatorError("evaluator label schema differs")
    if len(rows) != EXPECTED_PARTICIPANTS:
        raise NasEvaluatorError("evaluator participant count differs")
    row_ids = [row["row_id"] for row in rows]
    if len(set(row_ids)) != len(row_ids):
        raise NasEvaluatorError("evaluator row ID is duplicated")
    values = np.asarray([float(row["nash_crn_component_sum"]) for row in rows])
    if np.any(values < 0) or np.any(values > 8) or not np.all(np.isfinite(values)):
        raise NasEvaluatorError("an activity sum is off its 0-to-8 scale")
    return row_ids, values


def evaluate_gate(
    *, gate: Mapping[str, Any], candidate_id: str, comparator_id: str,
    candidate_metrics: Mapping[str, float], comparator_metrics: Mapping[str, float],
    paired_gain: Mapping[str, float], seed_summary: Mapping[str, Any],
    candidate_null: Mapping[str, float],
) -> dict[str, Any]:
    spec = gate["development_advance_conditions"]
    applicability = gate["seed_direction_applicability"]
    minimum_gain = float(spec["minimum_absolute_spearman_gain_over_strongest_baseline"])
    required_seeds = int(str(spec["minimum_seed_direction_consistency"]).split("_")[0])
    determinism_threshold = float(applicability["seed_determinism_correlation_threshold"])
    gain = float(candidate_metrics[PRIMARY_METRIC] - comparator_metrics[PRIMARY_METRIC])

    conditions: list[dict[str, Any]] = [
        {
            "condition_id": "minimum_absolute_spearman_gain_over_strongest_baseline",
            "requirement": f"Spearman gain >= {minimum_gain:.2f} over {comparator_id}",
            "observed": gain,
            "state": "met" if gain >= minimum_gain else "not_met",
        },
        {
            "condition_id": "gain_direction_required_bootstrap_lower_bound_above_zero",
            "requirement": "paired participant-bootstrap 2.5% bound on the Spearman gain > 0",
            "observed": float(paired_gain["ci_low"]),
            "state": "met" if paired_gain["ci_low"] > 0.0 else "not_met",
        },
    ]

    # A deterministic fit must not earn this condition.  Five identical vectors
    # are trivially all in the same direction, and a condition that pays out for
    # free inflates the headline without adding evidence.
    deterministic = bool(seed_summary["min_pairwise_correlation"] > determinism_threshold)
    seed_condition: dict[str, Any] = {
        "condition_id": "minimum_seed_direction_consistency",
        "requirement": f"{required_seeds} of {seed_summary['seeds']} seeds above {comparator_id}",
        "observed": int(seed_summary["seeds_above_comparator"]),
        "min_pairwise_correlation": float(seed_summary["min_pairwise_correlation"]),
        "determinism_threshold": determinism_threshold,
    }
    if deterministic:
        seed_condition["state"] = "not_applicable"
        seed_condition["reason"] = (
            f"the five seed prediction vectors correlate at a minimum of "
            f"{seed_summary['min_pairwise_correlation']:.6f}, above the {determinism_threshold} "
            "threshold, so this is one deterministic fit rather than five stability "
            "replicates; the condition is recorded not_applicable and is never met"
        )
    else:
        seed_condition["state"] = (
            "met" if seed_summary["seeds_above_comparator"] >= required_seeds else "not_met"
        )
    conditions.append(seed_condition)

    conditions.append({
        "condition_id": "candidate_must_exceed_its_own_permutation_null_p95",
        "requirement": "candidate Spearman above the 95th percentile of its own permutation null",
        "observed": float(candidate_metrics[PRIMARY_METRIC]),
        "null_p95": float(candidate_null["null_p95"]),
        "permutation_p_value": float(candidate_null["permutation_p_value"]),
        "state": (
            "met" if candidate_metrics[PRIMARY_METRIC] > candidate_null["null_p95"]
            else "not_met"
        ),
    })

    applicable = [c for c in conditions if c["state"] != "not_applicable"]
    met = [c for c in applicable if c["state"] == "met"]
    return {
        "gate_id": gate["gate_id"],
        "primary_endpoint": gate["primary_endpoint"],
        "candidate_model_id": candidate_id,
        "candidate_selection_basis": "source_only_pooled_oof_spearman_named_before_any_outcome_join",
        "comparator_model_id": comparator_id,
        "comparator_selection_basis": "strongest_non_candidate_baseline_on_the_target_post_hoc",
        "candidate_spearman": float(candidate_metrics[PRIMARY_METRIC]),
        "comparator_spearman": float(comparator_metrics[PRIMARY_METRIC]),
        "spearman_gain": gain,
        "paired_bootstrap_gain_ci_low": float(paired_gain["ci_low"]),
        "paired_bootstrap_gain_ci_high": float(paired_gain["ci_high"]),
        "conditions": conditions,
        "conditions_met": len(met),
        "conditions_applicable": len(applicable),
        "conditions_total": len(conditions),
        "conditions_not_applicable": [
            {"condition_id": c["condition_id"], "reason": c["reason"]}
            for c in conditions if c["state"] == "not_applicable"
        ],
        "verdict": "PASS" if len(met) == len(applicable) and applicable else "FAIL",
        "verdict_form": f"{len(met)} of {len(applicable)} applicable ({len(conditions)} total)",
        "non_gated_metrics_reported_only": list(NON_GATED_METRICS),
        "external_development_only": bool(gate["external_development_only"]),
        "champion_eligible": bool(gate["champion_eligible"]),
    }


def evaluate(
    *, labels_path: Path, prediction_roots: Sequence[Path], gate_path: Path,
    candidate_model_id: str, output: Path, bootstrap_replicates: int = 10000,
    bootstrap_seed: int = 20260826, permutation_replicates: int = 10000,
    permutation_seed: int = 20260826,
    prediction_bundle_artifacts_sha256: Sequence[str] | None = None,
) -> dict[str, Any]:
    if output.exists():
        raise NasEvaluatorError(f"refusing to overwrite evaluation: {output}")
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if (
        gate.get("primary_endpoint") != PRIMARY_METRIC
        or gate.get("external_development_only") is not True
        or gate.get("champion_eligible") is not False
        or "seed_direction_applicability" not in gate
    ):
        raise NasEvaluatorError("promotion gate differs")
    row_ids, observed = load_outcomes(labels_path)

    predictions: dict[str, np.ndarray] = {}
    seed_stacks: dict[str, np.ndarray] = {}
    seed_roster: dict[str, list[int]] = {}
    file_sha: dict[str, str] = {}
    artifacts_sha: dict[str, str] = {}
    for index, root in enumerate(prediction_roots):
        model_id, values = load_predictions(root / "predictions.tsv", row_ids)
        if model_id in predictions:
            raise NasEvaluatorError("model prediction bundle is duplicated")
        seeds, stack = load_seed_scores(root / "seed_scores.tsv", row_ids)
        if not np.allclose(stack.mean(axis=0), values, rtol=0.0, atol=1e-9):
            raise NasEvaluatorError("seed scores do not average to the frozen prediction")
        predictions[model_id] = values
        seed_stacks[model_id] = stack
        seed_roster[model_id] = seeds
        file_sha[model_id] = sha256_file(root / "predictions.tsv")
        if prediction_bundle_artifacts_sha256 is not None:
            artifacts_sha[model_id] = prediction_bundle_artifacts_sha256[index]
    missing = sorted(set(MANDATORY_BASELINES) - set(predictions))
    if missing:
        raise NasEvaluatorError(f"mandatory baselines are missing: {missing}")
    if candidate_model_id not in predictions:
        raise NasEvaluatorError("candidate model has no prediction bundle")

    output.mkdir(parents=True)
    point = {m: metrics(observed, v) for m, v in predictions.items()}
    nulls = {
        m: permutation_null(observed, v, replicates=permutation_replicates,
                            seed=permutation_seed)
        for m, v in sorted(predictions.items())
    }
    comparator_pool = [m for m in MANDATORY_BASELINES if m != candidate_model_id]
    comparator_id = min(
        comparator_pool,
        key=lambda v: (-point[v][PRIMARY_METRIC],
                       point[v]["participant_mae_nash_crn_component_sum"], v),
    )

    rng = np.random.default_rng(bootstrap_seed)
    n = len(observed)
    boot = {m: {k: np.empty(bootstrap_replicates) for k in METRIC_IDS} for m in predictions}
    gains = {m: np.empty(bootstrap_replicates) for m in predictions}
    usable = np.zeros(bootstrap_replicates, dtype=bool)
    for replicate in range(bootstrap_replicates):
        sampled = rng.integers(0, n, size=n)
        drawn = observed[sampled]
        if len(np.unique(drawn)) < 2:
            continue
        usable[replicate] = True
        comparator_value = metrics(drawn, predictions[comparator_id][sampled])
        for model_id, values in predictions.items():
            value = metrics(drawn, values[sampled])
            for key in METRIC_IDS:
                boot[model_id][key][replicate] = value[key]
            gains[model_id][replicate] = value[PRIMARY_METRIC] - comparator_value[PRIMARY_METRIC]
    if int(np.sum(usable)) < bootstrap_replicates // 2:
        raise NasEvaluatorError("participant bootstrap is degenerate")

    metric_rows, boot_rows = [], []
    for model_id in sorted(predictions):
        for key in METRIC_IDS:
            estimate = point[model_id][key]
            vector = boot[model_id][key][usable]
            vector = vector[np.isfinite(vector)]
            metric_rows.append({"model_id": model_id, "metric": key, "estimate": estimate,
                                "gated": key == PRIMARY_METRIC})
            boot_rows.append({
                "model_id": model_id, "metric": key, "estimate": estimate,
                "ci_low": float(np.quantile(vector, 0.025)) if vector.size else float("nan"),
                "ci_high": float(np.quantile(vector, 0.975)) if vector.size else float("nan"),
                "bootstrap_replicates_used": int(np.sum(usable)),
            })
    paired_rows, paired_by_model = [], {}
    for model_id in sorted(predictions):
        if model_id == comparator_id:
            continue
        vector = gains[model_id][usable]
        record = {
            "model_id": model_id, "comparator_model_id": comparator_id,
            "spearman_gain": float(point[model_id][PRIMARY_METRIC]
                                   - point[comparator_id][PRIMARY_METRIC]),
            "ci_low": float(np.quantile(vector, 0.025)),
            "ci_high": float(np.quantile(vector, 0.975)),
            "paired_draws_above_zero_fraction": float(np.mean(vector > 0.0)),
        }
        paired_rows.append(record)
        paired_by_model[model_id] = record

    stack = seed_stacks[candidate_model_id]
    seed_spearman = [spearman(observed, stack[i]) for i in range(stack.shape[0])]
    seed_summary = {
        "seeds": int(stack.shape[0]),
        "min_pairwise_correlation": min_pairwise_correlation(stack),
        "seed_spearman": [float(v) for v in seed_spearman],
        "max_pairwise_prediction_difference": float(
            np.max(stack.max(axis=0) - stack.min(axis=0))),
        "seeds_above_comparator": int(
            np.sum(np.asarray(seed_spearman) > point[comparator_id][PRIMARY_METRIC])),
    }
    verdict = evaluate_gate(
        gate=gate, candidate_id=candidate_model_id, comparator_id=comparator_id,
        candidate_metrics=point[candidate_model_id],
        comparator_metrics=point[comparator_id],
        paired_gain=paired_by_model[candidate_model_id],
        seed_summary=seed_summary, candidate_null=nulls[candidate_model_id],
    )

    null_rows = [{
        "model_id": m, "observed_spearman": point[m][PRIMARY_METRIC],
        **{k: nulls[m][k] for k in sorted(nulls[m])},
    } for m in sorted(predictions)]

    write_tsv(output / "metrics.tsv", ("model_id", "metric", "estimate", "gated"), metric_rows)
    write_tsv(output / "bootstrap_intervals.tsv",
              ("model_id", "metric", "estimate", "ci_low", "ci_high",
               "bootstrap_replicates_used"), boot_rows)
    write_tsv(output / "paired_vs_comparator.tsv", tuple(paired_rows[0]), paired_rows)
    write_tsv(output / "permutation_null.tsv", tuple(null_rows[0]), null_rows)
    (output / "promotion_gate_verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "seed_determinism.json").write_text(
        json.dumps(seed_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    receipt = {
        "schema_version": "masld-bench-gse267145-nas-sum-evaluation-v1",
        "status": "passed_external_development_scoring",
        "task_id": "gse267145_nas_sum_transfer",
        "participants": EXPECTED_PARTICIPANTS,
        "biological_unit": "participant",
        "target_activity_mean": float(np.mean(observed)),
        "target_activity_range": [float(np.min(observed)), float(np.max(observed))],
        "target_tied_zeros": int(np.sum(observed == 0)),
        "models": sorted(predictions),
        "mandatory_baselines": list(MANDATORY_BASELINES),
        "candidate_model_id": candidate_model_id,
        "comparator_model_id": comparator_id,
        "primary_metric": PRIMARY_METRIC,
        "non_gated_metrics": list(NON_GATED_METRICS),
        "non_gated_reason": "not invariant to the declared 1.688-point severity shift",
        "analytic_rank_null_sd": float(1.0 / np.sqrt(EXPECTED_PARTICIPANTS - 1)),
        "bootstrap_method": "unstratified_participant_resample_paired_across_models",
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_replicates_used": int(np.sum(usable)),
        "bootstrap_seed": bootstrap_seed,
        "permutation_replicates": permutation_replicates,
        "permutation_seed": permutation_seed,
        "promotion_gate_verdict": verdict["verdict"],
        "promotion_gate_verdict_form": verdict["verdict_form"],
        "conditions_met": verdict["conditions_met"],
        "conditions_applicable": verdict["conditions_applicable"],
        "conditions_total": verdict["conditions_total"],
        "endpoint_binarised": False,
        "recentred_or_rescaled_onto_the_target": False,
        "model_selected_from_target_outcomes": False,
        "prediction_file_sha256": dict(sorted(file_sha.items())),
        "prediction_bundle_artifacts_sha256": dict(sorted(artifacts_sha.items())),
        "external_development_only": True,
        "champion_claim_eligible": False,
        "absolute_calibration_claim_eligible": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluator-only", required=True, type=Path)
    parser.add_argument("--evaluator-only-artifacts-sha256", required=True)
    parser.add_argument("--gate", required=True, type=Path)
    parser.add_argument("--gate-sha256", required=True)
    parser.add_argument("--candidate-model-id", required=True)
    parser.add_argument("--prediction", action="append", required=True)
    parser.add_argument("--prediction-artifacts-sha256", action="append", required=True)
    parser.add_argument("--output", required=True, type=Path)
    a = parser.parse_args()
    if sha256_file(a.evaluator_only / "ARTIFACTS.json") != a.evaluator_only_artifacts_sha256:
        raise NasEvaluatorError("evaluator-only ARTIFACTS SHA differs")
    if sha256_file(a.gate) != a.gate_sha256:
        raise NasEvaluatorError("promotion gate SHA differs")
    roots = [Path(v) for v in a.prediction]
    for root, expected in zip(roots, a.prediction_artifacts_sha256, strict=True):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise NasEvaluatorError(f"prediction ARTIFACTS SHA differs: {root}")
    print(json.dumps(evaluate(
        labels_path=a.evaluator_only / "labels.tsv", prediction_roots=roots,
        gate_path=a.gate, candidate_model_id=a.candidate_model_id, output=a.output,
        prediction_bundle_artifacts_sha256=a.prediction_artifacts_sha256), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
