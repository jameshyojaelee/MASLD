#!/usr/bin/env python3
"""Independently join and evaluate frozen GSE267145 histology predictions."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np


STAGE3 = ("NOR", "NAFL", "NASH")
FIBROSIS_GROUP3 = ("F0", "F1", "F2_3")
MODEL_IDS = (
    "training_stage_distribution",
    "rna_hvg_pca_elastic_net",
    "rna_hvg_pca_linear_svm",
    "rna_hvg_pca_nearest_centroid",
    "rna_hvg_pca_knn",
    "h3_variance_pca_elastic_net",
    "h3_variance_pca_linear_svm",
    "h3_variance_pca_nearest_centroid",
    "h3_variance_pca_knn",
    "block_pca_elastic_net",
    "calibrated_late_fusion",
)
SEEDS = (1701, 1709, 1721, 1723, 1733)
PREDICTION_FIELDS = (
    "participant_id",
    "outer_fold",
    "probability_NOR",
    "probability_NAFL",
    "probability_NASH",
    "predicted_stage3",
    "predicted_nash_crn_component_sum",
    "predicted_fibrosis_cumulative_expected",
    "predicted_fibrosis_regression",
    "probability_fibrosis_F0",
    "probability_fibrosis_F1",
    "probability_fibrosis_F2_3",
    "predicted_fibrosis_group3",
)
METRIC_DIRECTIONS = {
    "stage3_macro_f1": "maximize",
    "stage3_multiclass_brier": "minimize",
    "fibrosis_group3_macro_f1": "maximize",
    "fibrosis_cumulative_ordinal_mae": "minimize",
    "fibrosis_cumulative_spearman": "maximize",
    "fibrosis_regression_mae": "minimize",
    "fibrosis_regression_spearman": "maximize",
    "nash_crn_component_sum_spearman": "maximize",
    "nash_crn_component_sum_mae": "minimize",
}
TASK_NATIVE_BASELINE_IDS = (
    "training_stage_distribution",
    "rna_hvg_pca_elastic_net",
    "rna_hvg_pca_linear_svm",
    "rna_hvg_pca_nearest_centroid",
    "rna_hvg_pca_knn",
    "h3_variance_pca_elastic_net",
    "h3_variance_pca_linear_svm",
    "h3_variance_pca_nearest_centroid",
    "h3_variance_pca_knn",
)


class HistologyEvaluationError(RuntimeError):
    """Raised when a frozen prediction/outcome join does not meet the requirements."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise HistologyEvaluationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _macro_f1(observed: np.ndarray, predicted: np.ndarray, roster: Sequence[str]) -> float:
    scores: list[float] = []
    for label in roster:
        tp = int(np.sum((observed == label) & (predicted == label)))
        fp = int(np.sum((observed != label) & (predicted == label)))
        fn = int(np.sum((observed == label) & (predicted != label)))
        denominator = 2 * tp + fp + fn
        scores.append(0.0 if denominator == 0 else 2.0 * tp / denominator)
    return float(np.mean(scores))


def _spearman(observed: np.ndarray, predicted: np.ndarray) -> float:
    from scipy.stats import spearmanr

    value = float(spearmanr(observed, predicted).statistic)
    return value


def _probabilities(rows: Sequence[Mapping[str, str]], prefix: str, roster: Sequence[str]) -> np.ndarray:
    values = np.asarray(
        [[float(row[f"{prefix}{label}"]) for label in roster] for row in rows],
        dtype=np.float64,
    )
    if (
        not np.all(np.isfinite(values))
        or np.any(values < 0)
        or not np.allclose(values.sum(axis=1), 1.0, atol=1e-6)
    ):
        raise HistologyEvaluationError("prediction probabilities are invalid")
    return values


def _brier(observed: np.ndarray, probabilities: np.ndarray, roster: Sequence[str]) -> float:
    one_hot = np.zeros_like(probabilities)
    lookup = {label: index for index, label in enumerate(roster)}
    for index, label in enumerate(observed):
        one_hot[index, lookup[str(label)]] = 1.0
    return float(np.mean(np.sum((probabilities - one_hot) ** 2, axis=1)))


def _fibrosis_group(value: int) -> str:
    return "F0" if value == 0 else "F1" if value == 1 else "F2_3"


def _metric_functions(
    outcomes: Mapping[str, np.ndarray], predictions: Mapping[str, np.ndarray]
) -> dict[str, Callable[[np.ndarray], float]]:
    return {
        "stage3_macro_f1": lambda index: _macro_f1(
            outcomes["stage3"][index], predictions["stage3"][index], STAGE3
        ),
        "stage3_multiclass_brier": lambda index: _brier(
            outcomes["stage3"][index], predictions["stage_probabilities"][index], STAGE3
        ),
        "fibrosis_group3_macro_f1": lambda index: _macro_f1(
            outcomes["fibrosis_group3"][index], predictions["fibrosis_group3"][index], FIBROSIS_GROUP3
        ),
        "fibrosis_cumulative_ordinal_mae": lambda index: float(
            np.mean(np.abs(outcomes["fibrosis"][index] - predictions["fibrosis_cumulative"][index]))
        ),
        "fibrosis_cumulative_spearman": lambda index: _spearman(
            outcomes["fibrosis"][index], predictions["fibrosis_cumulative"][index]
        ),
        "fibrosis_regression_mae": lambda index: float(
            np.mean(np.abs(outcomes["fibrosis"][index] - predictions["fibrosis_regression"][index]))
        ),
        "fibrosis_regression_spearman": lambda index: _spearman(
            outcomes["fibrosis"][index], predictions["fibrosis_regression"][index]
        ),
        "nash_crn_component_sum_spearman": lambda index: _spearman(
            outcomes["nas"][index], predictions["nas"][index]
        ),
        "nash_crn_component_sum_mae": lambda index: float(
            np.mean(np.abs(outcomes["nas"][index] - predictions["nas"][index]))
        ),
    }


def _bootstrap_interval(
    metric: Callable[[np.ndarray], float], bootstrap_indices: np.ndarray
) -> tuple[float | None, float | None, int]:
    values: list[float] = []
    for indices in bootstrap_indices:
        try:
            value = metric(indices)
        except (ValueError, KeyError, FloatingPointError):
            continue
        if math.isfinite(value):
            values.append(value)
    if not values:
        return None, None, 0
    return (
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
        len(values),
    )


def _paired_bootstrap_difference(
    model_metric: Callable[[np.ndarray], float],
    baseline_metric: Callable[[np.ndarray], float],
    bootstrap_indices: np.ndarray,
    *,
    direction: str,
) -> tuple[float | None, float | None, int, float | None]:
    sign = 1.0 if direction == "maximize" else -1.0
    values: list[float] = []
    for indices in bootstrap_indices:
        try:
            value = sign * (model_metric(indices) - baseline_metric(indices))
        except (ValueError, KeyError, FloatingPointError):
            continue
        if math.isfinite(value):
            values.append(value)
    if not values:
        return None, None, 0, None
    array = np.asarray(values, dtype=np.float64)
    return (
        float(np.quantile(array, 0.025)),
        float(np.quantile(array, 0.975)),
        len(values),
        float(np.mean(array > 0)),
    )


def _pareto_fronts(values: Mapping[str, Mapping[str, float]]) -> dict[str, int]:
    """Return nondominated front numbers with every endpoint oriented higher-is-better."""
    remaining = set(values)
    ranks: dict[str, int] = {}
    rank = 1
    while remaining:
        front: list[str] = []
        for candidate in sorted(remaining):
            dominated = False
            for other in remaining - {candidate}:
                other_values = values[other]
                candidate_values = values[candidate]
                if all(
                    other_values[metric] >= candidate_values[metric]
                    for metric in METRIC_DIRECTIONS
                ) and any(
                    other_values[metric] > candidate_values[metric]
                    for metric in METRIC_DIRECTIONS
                ):
                    dominated = True
                    break
            if not dominated:
                front.append(candidate)
        if not front:
            raise HistologyEvaluationError("Pareto ranking failed to identify a front")
        for model_id in front:
            ranks[model_id] = rank
            remaining.remove(model_id)
        rank += 1
    return ranks


def evaluate_campaign(
    *,
    predictions_root: Path,
    outcomes_path: Path,
    folds_path: Path,
    output: Path,
    bootstrap_replicates: int = 10_000,
    bootstrap_seed: int = 20260824,
) -> dict[str, Any]:
    if output.exists():
        raise HistologyEvaluationError(f"refusing to overwrite evaluation: {output}")
    endpoint_fields, endpoint_rows = read_tsv(outcomes_path)
    fold_fields, fold_rows = read_tsv(folds_path)
    required = {
        "participant_id",
        "outer_fold",
        "stage3",
        "stage5",
        "nash_crn_component_sum",
        "fibrosis",
        "recorded_sex",
    }
    if not required <= set(endpoint_fields) or fold_fields != ("participant_id", "outer_fold"):
        raise HistologyEvaluationError("endpoint or fold schema differs")
    participant_ids = [row["participant_id"] for row in endpoint_rows]
    if (
        len(participant_ids) != 99
        or participant_ids != [row["participant_id"] for row in fold_rows]
        or [row["outer_fold"] for row in endpoint_rows]
        != [row["outer_fold"] for row in fold_rows]
    ):
        raise HistologyEvaluationError("endpoint/fold participant axes differ")
    outcomes = {
        "stage3": np.asarray([row["stage3"] for row in endpoint_rows], dtype=str),
        "fibrosis": np.asarray([int(row["fibrosis"]) for row in endpoint_rows], dtype=float),
        "fibrosis_group3": np.asarray(
            [_fibrosis_group(int(row["fibrosis"])) for row in endpoint_rows], dtype=str
        ),
        "nas": np.asarray([int(row["nash_crn_component_sum"]) for row in endpoint_rows], dtype=float),
        "sex": np.asarray([row["recorded_sex"] for row in endpoint_rows], dtype=str),
    }
    if Counter(outcomes["sex"]) != {"F": 85, "M": 14} or set(
        outcomes["stage3"][outcomes["sex"] == "M"]
    ) != {"NASH"}:
        raise HistologyEvaluationError("recorded-sex limitation differs")
    if Counter(outcomes["fibrosis"])[3.0] != 4:
        raise HistologyEvaluationError("sparse fibrosis3 census differs")
    rng = np.random.default_rng(bootstrap_seed)
    bootstrap_indices = rng.integers(0, 99, size=(bootstrap_replicates, 99), endpoint=False)
    output.mkdir(parents=True)
    metric_rows: list[dict[str, Any]] = []
    sex_rows: list[dict[str, Any]] = []
    ensemble_functions: dict[str, dict[str, Callable[[np.ndarray], float]]] = {}
    ensemble_values: dict[str, dict[str, float]] = {}
    expected_files = [
        (model_id, None, predictions_root / f"{model_id}.tsv") for model_id in MODEL_IDS
    ] + [
        (model_id, seed, predictions_root / f"{model_id}--seed-{seed}.tsv")
        for model_id in MODEL_IDS
        for seed in SEEDS
    ]
    if set(predictions_root.glob("*.tsv")) != {path for _, _, path in expected_files}:
        raise HistologyEvaluationError("prediction file roster differs")
    for model_id, seed, path in expected_files:
        fields, rows = read_tsv(path)
        if fields != PREDICTION_FIELDS or len(rows) != 99:
            raise HistologyEvaluationError("prediction schema or row count differs")
        if (
            [row["participant_id"] for row in rows] != participant_ids
            or [row["outer_fold"] for row in rows]
            != [row["outer_fold"] for row in fold_rows]
        ):
            raise HistologyEvaluationError("prediction row/fold universe differs")
        stage_probabilities = _probabilities(rows, "probability_", STAGE3)
        group_probabilities = _probabilities(rows, "probability_fibrosis_", FIBROSIS_GROUP3)
        predictions = {
            "stage_probabilities": stage_probabilities,
            "stage3": np.asarray([row["predicted_stage3"] for row in rows], dtype=str),
            "fibrosis_group3": np.asarray(
                [row["predicted_fibrosis_group3"] for row in rows], dtype=str
            ),
            "fibrosis_cumulative": np.asarray(
                [float(row["predicted_fibrosis_cumulative_expected"]) for row in rows]
            ),
            "fibrosis_regression": np.asarray(
                [float(row["predicted_fibrosis_regression"]) for row in rows]
            ),
            "nas": np.asarray(
                [float(row["predicted_nash_crn_component_sum"]) for row in rows]
            ),
        }
        if (
            not set(predictions["stage3"]) <= set(STAGE3)
            or not set(predictions["fibrosis_group3"]) <= set(FIBROSIS_GROUP3)
            or any(not np.all(np.isfinite(predictions[key])) for key in ("fibrosis_cumulative", "fibrosis_regression", "nas"))
        ):
            raise HistologyEvaluationError("prediction values differ")
        functions = _metric_functions(outcomes, predictions)
        full_index = np.arange(99)
        if seed is None:
            ensemble_functions[model_id] = functions
            ensemble_values[model_id] = {
                metric_id: function(full_index)
                for metric_id, function in functions.items()
            }
        for metric_id, function in functions.items():
            value = function(full_index)
            low, high, valid = _bootstrap_interval(function, bootstrap_indices)
            metric_rows.append(
                {
                    "model_id": model_id,
                    "seed": "ensemble" if seed is None else seed,
                    "metric_id": metric_id,
                    "value": format(value, ".17g") if math.isfinite(value) else "not_estimable",
                    "ci95_low": "not_estimable" if low is None else format(low, ".17g"),
                    "ci95_high": "not_estimable" if high is None else format(high, ".17g"),
                    "valid_bootstrap_replicates": valid,
                    "biological_resampling_unit": "participant",
                }
            )
        for recorded_sex in ("F", "M"):
            mask = outcomes["sex"] == recorded_sex
            nash_mask = mask & (outcomes["stage3"] == "NASH")
            sex_rows.append(
                {
                    "model_id": model_id,
                    "seed": "ensemble" if seed is None else seed,
                    "recorded_sex": recorded_sex,
                    "participants": int(mask.sum()),
                    "stage3_census": canonical_json(dict(sorted(Counter(outcomes["stage3"][mask]).items()))),
                    "stage3_brier": format(
                        _brier(outcomes["stage3"][mask], stage_probabilities[mask], STAGE3), ".17g"
                    ),
                    "stage3_macro_f1": (
                        format(_macro_f1(outcomes["stage3"][mask], predictions["stage3"][mask], STAGE3), ".17g")
                        if recorded_sex == "F"
                        else "not_applicable_missing_NOR_and_NAFL"
                    ),
                    "nash_recall": format(
                        float(np.mean(predictions["stage3"][nash_mask] == "NASH")), ".17g"
                    ),
                    "mean_nash_true_class_probability": format(
                        float(np.mean(stage_probabilities[nash_mask, STAGE3.index("NASH")])), ".17g"
                    ),
                    "between_sex_global_gap_claim_allowed": "false",
                }
            )
    write_tsv(
        output / "metrics.tsv",
        (
            "model_id",
            "seed",
            "metric_id",
            "value",
            "ci95_low",
            "ci95_high",
            "valid_bootstrap_replicates",
            "biological_resampling_unit",
        ),
        metric_rows,
    )
    write_tsv(
        output / "sex_error_audit.tsv",
        (
            "model_id",
            "seed",
            "recorded_sex",
            "participants",
            "stage3_census",
            "stage3_brier",
            "stage3_macro_f1",
            "nash_recall",
            "mean_nash_true_class_probability",
            "between_sex_global_gap_claim_allowed",
        ),
        sex_rows,
    )
    strongest_by_metric: dict[str, str] = {}
    for metric_id, direction in METRIC_DIRECTIONS.items():
        oriented = [
            (
                ensemble_values[model_id][metric_id]
                if direction == "maximize"
                else -ensemble_values[model_id][metric_id],
                model_id,
            )
            for model_id in TASK_NATIVE_BASELINE_IDS
        ]
        strongest_by_metric[metric_id] = min(
            oriented, key=lambda item: (-item[0], item[1])
        )[1]
    difference_rows: list[dict[str, Any]] = []
    oriented_values: dict[str, dict[str, float]] = {}
    for model_id in MODEL_IDS:
        oriented_values[model_id] = {}
        for metric_id, direction in METRIC_DIRECTIONS.items():
            baseline_id = strongest_by_metric[metric_id]
            model_value = ensemble_values[model_id][metric_id]
            baseline_value = ensemble_values[baseline_id][metric_id]
            sign = 1.0 if direction == "maximize" else -1.0
            difference = sign * (model_value - baseline_value)
            low, high, valid, positive_fraction = _paired_bootstrap_difference(
                ensemble_functions[model_id][metric_id],
                ensemble_functions[baseline_id][metric_id],
                bootstrap_indices,
                direction=direction,
            )
            oriented_values[model_id][metric_id] = sign * model_value
            difference_rows.append(
                {
                    "model_id": model_id,
                    "metric_id": metric_id,
                    "direction": direction,
                    "strongest_task_native_baseline_id": baseline_id,
                    "model_value": format(model_value, ".17g"),
                    "baseline_value": format(baseline_value, ".17g"),
                    "oriented_difference_positive_is_better": format(difference, ".17g"),
                    "paired_ci95_low": "not_estimable" if low is None else format(low, ".17g"),
                    "paired_ci95_high": "not_estimable" if high is None else format(high, ".17g"),
                    "valid_bootstrap_replicates": valid,
                    "bootstrap_fraction_strictly_positive": (
                        "not_estimable"
                        if positive_fraction is None
                        else format(positive_fraction, ".17g")
                    ),
                    "confirmatory_inference_allowed": "false",
                }
            )
    write_tsv(
        output / "paired_differences.tsv",
        (
            "model_id",
            "metric_id",
            "direction",
            "strongest_task_native_baseline_id",
            "model_value",
            "baseline_value",
            "oriented_difference_positive_is_better",
            "paired_ci95_low",
            "paired_ci95_high",
            "valid_bootstrap_replicates",
            "bootstrap_fraction_strictly_positive",
            "confirmatory_inference_allowed",
        ),
        difference_rows,
    )
    pareto_ranks = _pareto_fronts(oriented_values)
    pareto_rows: list[dict[str, Any]] = []
    for model_id in MODEL_IDS:
        favorable = sum(
            oriented_values[model_id][metric_id]
            >= oriented_values[strongest_by_metric[metric_id]][metric_id]
            for metric_id in METRIC_DIRECTIONS
        )
        pareto_rows.append(
            {
                "model_id": model_id,
                "pareto_front": pareto_ranks[model_id],
                "endpoints_at_or_above_strongest_task_native_baseline": favorable,
                "endpoint_count": len(METRIC_DIRECTIONS),
                "ranking_scope": "ensemble_predictions_single_cohort_development_only",
            }
        )
    write_tsv(
        output / "pareto_ranking.tsv",
        (
            "model_id",
            "pareto_front",
            "endpoints_at_or_above_strongest_task_native_baseline",
            "endpoint_count",
            "ranking_scope",
        ),
        sorted(pareto_rows, key=lambda row: (row["pareto_front"], row["model_id"])),
    )
    disposition = {
        "schema_version": "masld-bench-gse267145-histology-evaluator-disposition-v1",
        "status": "single_cohort_development_only",
        "primary_endpoint": "stage3_macro_f1",
        "fibrosis": {
            "limitation_code": "sparse_f3_inner_training_instability",
            "recorded_f3_participants": 4,
            "minimum_inner_training_f3": 1,
            "cumulative_ge3_instability_expected": True,
            "f3_specific_superiority_claim_allowed": False,
            "champion_gate_eligible": False,
            "exact_regression_and_cumulative_heads_reported_separately": True,
            "role": "secondary_instability_report_required",
        },
        "recorded_sex": {
            "limitation_code": "recorded_male_perfectly_confound_with_nash",
            "female_participants": 85,
            "male_participants": 14,
            "male_stage3_classes": ["NASH"],
            "male_stage3_macro_f1": "not_applicable_missing_NOR_and_NAFL",
            "global_between_sex_performance_gap_claim_allowed": False,
            "sex_fairness_claim_allowed": False,
            "champion_gate_eligible": False,
            "role": "descriptive_error_audit_only",
        },
        "participant_age": "structurally_missing_no_audit",
        "source_stage5": "descriptive_only_not_predicted_or_used_for_selection",
        "seeds_are_biological_replicates": False,
        "secondary_p_values_calculated": False,
        "bh_multiplicity_activated": False,
        "clinical_external_or_champion_claim_allowed": False,
        "strongest_baseline_selection_scope": "same_outer_predictions_descriptive_development_only",
        "paired_bootstrap_differences_confirmatory": False,
        "pareto_ranking_confirmatory": False,
    }
    (output / "disposition.json").write_text(
        json.dumps(disposition, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema_version": "masld-bench-gse267145-histology-independent-evaluation-v1",
        "status": "passed_development_metrics",
        "evaluator_id": "gse267145_participant_histology_v1",
        "models": len(MODEL_IDS),
        "seed_prediction_sets": len(MODEL_IDS) * len(SEEDS),
        "ensemble_prediction_sets": len(MODEL_IDS),
        "metrics_per_prediction_set": len(METRIC_DIRECTIONS),
        "paired_model_difference_rows": len(MODEL_IDS) * len(METRIC_DIRECTIONS),
        "pareto_ranked_models": len(MODEL_IDS),
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_unit": "participant",
        "prediction_bundle_mutated": False,
        "fit_source_imported": False,
        "source_stage5_used_for_selection": False,
        "recorded_sex_used_for_selection": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--predictions-artifacts-sha256", required=True)
    parser.add_argument("--outcomes", required=True, type=Path)
    parser.add_argument("--outcomes-artifacts-sha256", required=True)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--folds-artifacts-sha256", required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    for root, expected in (
        (arguments.predictions, arguments.predictions_artifacts_sha256),
        (arguments.outcomes, arguments.outcomes_artifacts_sha256),
        (arguments.folds, arguments.folds_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise HistologyEvaluationError(f"input ARTIFACTS SHA differs: {root}")
    result = evaluate_campaign(
        predictions_root=arguments.predictions,
        outcomes_path=arguments.outcomes / "participant_endpoints.tsv",
        folds_path=arguments.folds / "participant_outer_folds.tsv",
        output=arguments.output,
        bootstrap_replicates=arguments.bootstrap_replicates,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
