#!/usr/bin/env python3
"""Independently score frozen GSE49541 fibrosis-group transfer predictions.

Four things this evaluator refuses to do the easy way:

1.  Average precision is not scored against prevalence.  A random scorer on 72
    participants with 32 positives does not average 0.444, and its 95th
    percentile sits well above that, so a prevalence line manufactures
    significance.  The null is built by permuting the outcome against each
    model's own realised score vector, because a tied or low-cardinality scorer
    has a tighter null than a continuous one.  Both references are reported.
2.  Intervals come from a participant bootstrap, and the candidate-versus-
    baseline comparison is paired on the same resampled participants.
    Overlapping marginal intervals are not a test of no difference.
3.  Seed direction consistency is computed from the seed-level scores that the
    predictor actually emitted.  Where five schema seeds produced one distinct
    prediction vector, that is reported as one deterministic fit, not as five
    consistent replicates.
4.  The strongest comparator is chosen among the baselines other than the
    candidate, and the candidate was named by source-only evidence before any
    outcome was joined.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


CLASSES = ("mild_f0_f1", "advanced_f3_f4")
POSITIVE_CLASS = "advanced_f3_f4"
EXPECTED_GROUP_COUNTS = {"mild_f0_f1": 40, "advanced_f3_f4": 32}
EXPECTED_PARTICIPANTS = 72
MANDATORY_BASELINES = (
    "gene_median_elastic_net",
    "gene_median_linear_svm",
    "gene_median_pca_elastic_net",
    "per_array_rank_elastic_net",
    "training_prevalence",
)
METRIC_IDS = (
    "participant_auprc_advanced_f3_f4",
    "participant_auroc_advanced_f3_f4",
    "participant_macro_f1_fibrosis_group",
    "participant_brier_fibrosis_group",
)


class ExternalTransferEvaluatorError(RuntimeError):
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
            raise ExternalTransferEvaluatorError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
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


# ---------------------------------------------------------------------------
# Metrics.
# ---------------------------------------------------------------------------


def average_precision(observed: np.ndarray, scores: np.ndarray) -> float:
    from sklearn.metrics import average_precision_score

    return float(average_precision_score(observed, scores))


def area_under_roc(observed: np.ndarray, scores: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score

    if len(np.unique(observed)) != 2:
        return float("nan")
    return float(roc_auc_score(observed, scores))


def macro_f1(observed: np.ndarray, scores: np.ndarray, *, threshold: float = 0.5) -> float:
    predicted = (np.asarray(scores, dtype=np.float64) > threshold).astype(np.int64)
    per_class = []
    for class_id in (0, 1):
        true = observed == class_id
        called = predicted == class_id
        tp = int(np.sum(true & called))
        fp = int(np.sum(~true & called))
        fn = int(np.sum(true & ~called))
        per_class.append(0.0 if 2 * tp + fp + fn == 0 else 2 * tp / (2 * tp + fp + fn))
    return float(np.mean(per_class))


def brier(observed: np.ndarray, scores: np.ndarray) -> float:
    """Binary Brier score on the advanced-fibrosis indicator (lower is better)."""

    return float(np.mean((np.asarray(scores, dtype=np.float64) - observed) ** 2))


def metrics(observed: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    return {
        "participant_auprc_advanced_f3_f4": average_precision(observed, scores),
        "participant_auroc_advanced_f3_f4": area_under_roc(observed, scores),
        "participant_macro_f1_fibrosis_group": macro_f1(observed, scores),
        "participant_brier_fibrosis_group": brier(observed, scores),
    }


def class_statistics(observed: np.ndarray, scores: np.ndarray) -> list[dict[str, Any]]:
    predicted = (np.asarray(scores, dtype=np.float64) > 0.5).astype(np.int64)
    rows: list[dict[str, Any]] = []
    for class_id, label in enumerate(CLASSES):
        true = observed == class_id
        called = predicted == class_id
        tp = int(np.sum(true & called))
        fp = int(np.sum(~true & called))
        fn = int(np.sum(true & ~called))
        rows.append(
            {
                "class": label,
                "support": int(np.sum(true)),
                "precision": 0.0 if tp + fp == 0 else tp / (tp + fp),
                "recall": 0.0 if tp + fn == 0 else tp / (tp + fn),
                "f1": 0.0 if 2 * tp + fp + fn == 0 else 2 * tp / (2 * tp + fp + fn),
                "tp": tp,
                "fp": fp,
                "fn": fn,
            }
        )
    return rows


def permutation_null(
    observed: np.ndarray,
    scores: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, float]:
    """Null for a fixed score vector, built by permuting the outcome.

    Holding the scores fixed keeps the null's tie structure identical to the
    scorer's own, which is what makes a constant scorer's null narrow and a
    continuous scorer's null wide.  Comparing a continuous model against a tied
    scorer's null is how a false positive gets manufactured.
    """

    rng = np.random.default_rng(seed)
    labels = np.asarray(observed, dtype=np.int64).copy()
    values = np.empty(replicates, dtype=np.float64)
    roc = np.empty(replicates, dtype=np.float64)
    for index in range(replicates):
        rng.shuffle(labels)
        values[index] = average_precision(labels, scores)
        roc[index] = area_under_roc(labels, scores)
    point = average_precision(observed, scores)
    return {
        "null_replicates": int(replicates),
        "auprc_null_mean": float(np.mean(values)),
        "auprc_null_p50": float(np.quantile(values, 0.50)),
        "auprc_null_p95": float(np.quantile(values, 0.95)),
        "auprc_null_p99": float(np.quantile(values, 0.99)),
        "auprc_permutation_p_value": float(
            (1.0 + float(np.sum(values >= point))) / (replicates + 1.0)
        ),
        "auroc_null_mean": float(np.mean(roc)),
        "auroc_null_p95": float(np.quantile(roc, 0.95)),
        "distinct_score_values": int(len(np.unique(np.round(scores, 12)))),
    }


# ---------------------------------------------------------------------------
# Prediction and outcome loading.
# ---------------------------------------------------------------------------


def load_predictions(path: Path, row_ids: Sequence[str]) -> tuple[str, np.ndarray]:
    fields, rows = read_tsv(path)
    required = {
        "row_id",
        "model_id",
        "predicted_fibrosis_group",
        "probability_mild_f0_f1",
        "probability_advanced_f3_f4",
    }
    if not required <= set(fields) or len(rows) != len(row_ids):
        raise ExternalTransferEvaluatorError("prediction schema or row count differs")
    model_ids = {row["model_id"] for row in rows}
    if len(model_ids) != 1:
        raise ExternalTransferEvaluatorError("prediction bundle contains multiple models")
    by_id = {row["row_id"]: row for row in rows}
    if len(by_id) != len(rows) or set(by_id) != set(row_ids):
        raise ExternalTransferEvaluatorError("prediction row IDs differ")
    ordered = [by_id[value] for value in row_ids]
    positive = np.asarray(
        [float(row["probability_advanced_f3_f4"]) for row in ordered], dtype=np.float64
    )
    negative = np.asarray(
        [float(row["probability_mild_f0_f1"]) for row in ordered], dtype=np.float64
    )
    if (
        not np.all(np.isfinite(positive))
        or np.any(positive < 0.0)
        or np.any(positive > 1.0)
        or not np.allclose(positive + negative, 1.0, rtol=0.0, atol=1e-8)
    ):
        raise ExternalTransferEvaluatorError("prediction probabilities are invalid")
    expected = [CLASSES[1] if value > 0.5 else CLASSES[0] for value in positive]
    if expected != [row["predicted_fibrosis_group"] for row in ordered]:
        raise ExternalTransferEvaluatorError(
            "predicted labels differ from the fixed-threshold call"
        )
    return next(iter(model_ids)), positive


def load_seed_scores(
    path: Path, row_ids: Sequence[str]
) -> tuple[list[int], np.ndarray]:
    fields, rows = read_tsv(path)
    required = {"row_id", "model_seed", "probability_advanced_f3_f4"}
    if not required <= set(fields):
        raise ExternalTransferEvaluatorError("seed-score schema differs")
    seeds = sorted({int(row["model_seed"]) for row in rows})
    table = {
        (row["row_id"], int(row["model_seed"])): float(
            row["probability_advanced_f3_f4"]
        )
        for row in rows
    }
    if len(table) != len(rows) or len(rows) != len(seeds) * len(row_ids):
        raise ExternalTransferEvaluatorError("seed-score table is not rectangular")
    matrix = np.asarray(
        [[table[(row_id, seed)] for row_id in row_ids] for seed in seeds],
        dtype=np.float64,
    )
    if not np.all(np.isfinite(matrix)):
        raise ExternalTransferEvaluatorError("seed scores are not finite")
    return seeds, matrix


def load_outcomes(path: Path) -> tuple[list[str], np.ndarray]:
    fields, rows = read_tsv(path)
    if not {"row_id", "fibrosis_stage_group"} <= set(fields):
        raise ExternalTransferEvaluatorError("evaluator label schema differs")
    if len(rows) != EXPECTED_PARTICIPANTS:
        raise ExternalTransferEvaluatorError("evaluator participant count differs")
    row_ids = [row["row_id"] for row in rows]
    if len(set(row_ids)) != len(row_ids):
        raise ExternalTransferEvaluatorError("evaluator row ID is duplicated")
    groups = [row["fibrosis_stage_group"] for row in rows]
    if any(value not in CLASSES for value in groups):
        raise ExternalTransferEvaluatorError("evaluator fibrosis group differs")
    if Counter(groups) != Counter(EXPECTED_GROUP_COUNTS):
        raise ExternalTransferEvaluatorError(
            "external fibrosis-group counts differ from 40 mild and 32 advanced"
        )
    observed = np.asarray([CLASSES.index(value) for value in groups], dtype=np.int64)
    return row_ids, observed


# ---------------------------------------------------------------------------
# Check.
# ---------------------------------------------------------------------------


def seed_direction_consistency(
    observed: np.ndarray,
    seed_scores: np.ndarray,
    baseline_auprc: float,
) -> dict[str, Any]:
    """How many seeds put the candidate above the comparator, honestly counted."""

    distinct = int(len(np.unique(np.round(seed_scores, 12), axis=0)))
    per_seed = [
        average_precision(observed, seed_scores[index])
        for index in range(seed_scores.shape[0])
    ]
    above = int(np.sum(np.asarray(per_seed) > baseline_auprc))
    return {
        "seeds": int(seed_scores.shape[0]),
        "distinct_seed_prediction_vectors": distinct,
        "seed_auprc": [float(value) for value in per_seed],
        "seeds_above_comparator": above,
        "deterministic_single_fit": distinct == 1,
        "consistency_is_evidence_of_stability": distinct > 1,
        "interpretation": (
            "five schema seeds produced one identical fitted prediction, so the "
            "count is 5-of-5 by construction and carries no stability evidence"
            if distinct == 1
            else f"{distinct} distinct seed prediction vectors; the count is a real "
            "stability replicate"
        ),
    }


def evaluate_gate(
    *,
    gate: Mapping[str, Any],
    candidate_id: str,
    comparator_id: str,
    candidate_metrics: Mapping[str, float],
    comparator_metrics: Mapping[str, float],
    paired_gain: Mapping[str, float],
    seed_summary: Mapping[str, Any],
) -> dict[str, Any]:
    conditions_spec = gate["development_advance_conditions"]
    minimum_gain = float(conditions_spec["minimum_absolute_auprc_gain_over_strongest_baseline"])
    maximum_brier = float(conditions_spec["brier_score_degradation_maximum"])
    required_seeds = int(str(conditions_spec["minimum_seed_direction_consistency"]).split("_")[0])
    gain = float(
        candidate_metrics["participant_auprc_advanced_f3_f4"]
        - comparator_metrics["participant_auprc_advanced_f3_f4"]
    )
    brier_change = float(
        candidate_metrics["participant_brier_fibrosis_group"]
        - comparator_metrics["participant_brier_fibrosis_group"]
    )
    conditions = [
        {
            "condition_id": "minimum_absolute_auprc_gain_over_strongest_baseline",
            "requirement": f"AUPRC gain >= {minimum_gain:.2f} over {comparator_id}",
            "observed": gain,
            "met": bool(gain >= minimum_gain),
        },
        {
            "condition_id": "gain_direction_required_bootstrap_lower_bound_above_zero",
            "requirement": "paired participant-bootstrap 2.5% bound on the AUPRC gain > 0",
            "observed": float(paired_gain["ci_low"]),
            "met": bool(paired_gain["ci_low"] > 0.0),
        },
        {
            "condition_id": "minimum_seed_direction_consistency",
            "requirement": f"{required_seeds} of {seed_summary['seeds']} seeds above {comparator_id}",
            "observed": int(seed_summary["seeds_above_comparator"]),
            "met": bool(
                seed_summary["seeds_above_comparator"] >= required_seeds
                and seed_summary["consistency_is_evidence_of_stability"]
            ),
            "note": seed_summary["interpretation"],
        },
        {
            "condition_id": "brier_score_degradation_maximum",
            "requirement": f"Brier increase over {comparator_id} <= {maximum_brier:.2f}",
            "observed": brier_change,
            "met": bool(brier_change <= maximum_brier),
        },
    ]
    passed = all(condition["met"] for condition in conditions)
    return {
        "gate_id": gate["gate_id"],
        "primary_endpoint": gate["primary_endpoint"],
        "candidate_model_id": candidate_id,
        "candidate_selection_basis": "source_only_pooled_oof_average_precision_named_before_label_join",
        "comparator_model_id": comparator_id,
        "comparator_selection_basis": "strongest_non_candidate_baseline_on_the_external_outcome_post_hoc",
        "candidate_auprc": float(candidate_metrics["participant_auprc_advanced_f3_f4"]),
        "comparator_auprc": float(comparator_metrics["participant_auprc_advanced_f3_f4"]),
        "auprc_gain": gain,
        "paired_bootstrap_gain_ci_low": float(paired_gain["ci_low"]),
        "paired_bootstrap_gain_ci_high": float(paired_gain["ci_high"]),
        "brier_change": brier_change,
        "conditions": conditions,
        "conditions_met": int(sum(1 for condition in conditions if condition["met"])),
        "conditions_total": len(conditions),
        "verdict": "PASS" if passed else "FAIL",
        "external_development_only": bool(gate["external_development_only"]),
        "champion_eligible": bool(gate["champion_eligible"]),
    }


def evaluate(
    *,
    labels_path: Path,
    prediction_roots: Sequence[Path],
    gate_path: Path,
    candidate_model_id: str,
    output: Path,
    arm_id: str,
    gate_eligible_arm: bool,
    bootstrap_replicates: int = 10000,
    bootstrap_seed: int = 20260826,
    permutation_replicates: int = 10000,
    permutation_seed: int = 20260826,
    unfittable_baselines: Sequence[str] | None = None,
    prediction_bundle_artifacts_sha256: Sequence[str] | None = None,
) -> dict[str, Any]:
    if output.exists():
        raise ExternalTransferEvaluatorError(f"refusing to overwrite evaluation: {output}")
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if (
        gate.get("primary_endpoint") != "participant_auprc_advanced_f3_f4"
        or gate.get("external_development_only") is not True
        or gate.get("champion_eligible") is not False
    ):
        raise ExternalTransferEvaluatorError("promotion gate differs")
    row_ids, observed = load_outcomes(labels_path)
    prevalence = float(np.mean(observed))

    model_scores: dict[str, np.ndarray] = {}
    model_seed_scores: dict[str, np.ndarray] = {}
    model_seed_roster: dict[str, list[int]] = {}
    prediction_file_sha256: dict[str, str] = {}
    prediction_artifacts_sha256: dict[str, str] = {}
    if prediction_bundle_artifacts_sha256 is not None and len(
        prediction_bundle_artifacts_sha256
    ) != len(prediction_roots):
        raise ExternalTransferEvaluatorError("prediction artifact SHA roster differs")
    for index, root in enumerate(prediction_roots):
        model_id, scores = load_predictions(root / "predictions.tsv", row_ids)
        if model_id in model_scores:
            raise ExternalTransferEvaluatorError("model prediction bundle is duplicated")
        model_scores[model_id] = scores
        seeds, seed_matrix = load_seed_scores(root / "seed_scores.tsv", row_ids)
        if not np.allclose(seed_matrix.mean(axis=0), scores, rtol=0.0, atol=1e-9):
            raise ExternalTransferEvaluatorError(
                "seed scores do not average to the frozen prediction"
            )
        model_seed_scores[model_id] = seed_matrix
        model_seed_roster[model_id] = seeds
        prediction_file_sha256[model_id] = sha256_file(root / "predictions.tsv")
        if prediction_bundle_artifacts_sha256 is not None:
            prediction_artifacts_sha256[model_id] = prediction_bundle_artifacts_sha256[index]
    missing = sorted(set(MANDATORY_BASELINES) - set(model_scores))
    # A mandatory baseline may be absent only where the source fit recorded it as
    # unfittable and the caller declared it here.  Anything else is a dropped
    # comparator, and dropping a comparator is how a weak result gets to look
    # strong.
    declared = sorted(unfittable_baselines or ())
    if missing != declared:
        raise ExternalTransferEvaluatorError(
            f"mandatory baselines are missing without a recorded unfittable "
            f"declaration: missing={missing} declared_unfittable={declared}"
        )
    if declared and candidate_model_id in declared:
        raise ExternalTransferEvaluatorError("the candidate model is declared unfittable")
    if candidate_model_id not in model_scores:
        raise ExternalTransferEvaluatorError("candidate model has no prediction bundle")

    output.mkdir(parents=True)
    point_metrics = {
        model_id: metrics(observed, scores) for model_id, scores in model_scores.items()
    }
    null_rows: list[dict[str, Any]] = []
    null_by_model: dict[str, dict[str, float]] = {}
    for model_id in sorted(model_scores):
        null = permutation_null(
            observed,
            model_scores[model_id],
            replicates=permutation_replicates,
            seed=permutation_seed,
        )
        null_by_model[model_id] = null
        null_rows.append(
            {
                "model_id": model_id,
                "observed_auprc": point_metrics[model_id]["participant_auprc_advanced_f3_f4"],
                "prevalence_reference_not_a_null": prevalence,
                **{key: null[key] for key in sorted(null)},
            }
        )

    comparator_pool = [
        model_id
        for model_id in MANDATORY_BASELINES
        if model_id != candidate_model_id and model_id in model_scores
    ]
    if not comparator_pool:
        raise ExternalTransferEvaluatorError("no comparator baseline survives")
    comparator_id = min(
        comparator_pool,
        key=lambda value: (
            -point_metrics[value]["participant_auprc_advanced_f3_f4"],
            point_metrics[value]["participant_brier_fibrosis_group"],
            value,
        ),
    )

    rng = np.random.default_rng(bootstrap_seed)
    n = len(observed)
    bootstrap_values = {
        model_id: {metric_id: np.empty(bootstrap_replicates) for metric_id in METRIC_IDS}
        for model_id in model_scores
    }
    gain_draws = {
        model_id: np.empty(bootstrap_replicates) for model_id in model_scores
    }
    brier_draws = {
        model_id: np.empty(bootstrap_replicates) for model_id in model_scores
    }
    usable = np.zeros(bootstrap_replicates, dtype=bool)
    for replicate in range(bootstrap_replicates):
        sampled = rng.integers(0, n, size=n)
        drawn = observed[sampled]
        if len(np.unique(drawn)) != 2:
            continue
        usable[replicate] = True
        comparator_value = metrics(drawn, model_scores[comparator_id][sampled])
        for model_id, scores in model_scores.items():
            value = metrics(drawn, scores[sampled])
            for metric_id in METRIC_IDS:
                bootstrap_values[model_id][metric_id][replicate] = value[metric_id]
            gain_draws[model_id][replicate] = (
                value["participant_auprc_advanced_f3_f4"]
                - comparator_value["participant_auprc_advanced_f3_f4"]
            )
            brier_draws[model_id][replicate] = (
                value["participant_brier_fibrosis_group"]
                - comparator_value["participant_brier_fibrosis_group"]
            )
    if int(np.sum(usable)) < bootstrap_replicates // 2:
        raise ExternalTransferEvaluatorError("participant bootstrap is degenerate")

    metric_rows: list[dict[str, Any]] = []
    bootstrap_rows: list[dict[str, Any]] = []
    for model_id in sorted(model_scores):
        for metric_id in METRIC_IDS:
            estimate = point_metrics[model_id][metric_id]
            vector = bootstrap_values[model_id][metric_id][usable]
            metric_rows.append(
                {"model_id": model_id, "metric": metric_id, "estimate": estimate}
            )
            bootstrap_rows.append(
                {
                    "model_id": model_id,
                    "metric": metric_id,
                    "estimate": estimate,
                    "ci_low": float(np.quantile(vector, 0.025)),
                    "ci_high": float(np.quantile(vector, 0.975)),
                    "bootstrap_replicates_used": int(np.sum(usable)),
                }
            )
    paired_rows: list[dict[str, Any]] = []
    paired_by_model: dict[str, dict[str, float]] = {}
    for model_id in sorted(model_scores):
        if model_id == comparator_id:
            continue
        gains = gain_draws[model_id][usable]
        briers = brier_draws[model_id][usable]
        record = {
            "model_id": model_id,
            "comparator_model_id": comparator_id,
            "auprc_gain": float(
                point_metrics[model_id]["participant_auprc_advanced_f3_f4"]
                - point_metrics[comparator_id]["participant_auprc_advanced_f3_f4"]
            ),
            "ci_low": float(np.quantile(gains, 0.025)),
            "ci_high": float(np.quantile(gains, 0.975)),
            "paired_draws_above_zero_fraction": float(np.mean(gains > 0.0)),
            "brier_change": float(
                point_metrics[model_id]["participant_brier_fibrosis_group"]
                - point_metrics[comparator_id]["participant_brier_fibrosis_group"]
            ),
            "brier_change_ci_low": float(np.quantile(briers, 0.025)),
            "brier_change_ci_high": float(np.quantile(briers, 0.975)),
        }
        paired_rows.append(record)
        paired_by_model[model_id] = record

    seed_summary = seed_direction_consistency(
        observed,
        model_seed_scores[candidate_model_id],
        point_metrics[comparator_id]["participant_auprc_advanced_f3_f4"],
    )
    verdict = evaluate_gate(
        gate=gate,
        candidate_id=candidate_model_id,
        comparator_id=comparator_id,
        candidate_metrics=point_metrics[candidate_model_id],
        comparator_metrics=point_metrics[comparator_id],
        paired_gain=paired_by_model[candidate_model_id],
        seed_summary=seed_summary,
    )
    verdict["arm_id"] = arm_id
    verdict["gate_eligible_arm"] = bool(gate_eligible_arm)
    verdict["unfittable_baselines_declared"] = declared
    verdict["scored_baselines"] = sorted(set(MANDATORY_BASELINES) & set(model_scores))
    if not gate_eligible_arm:
        verdict["verdict_binding"] = False
        verdict["verdict_note"] = (
            "sensitivity arm: recorded for interpretation only and cannot advance the lane"
        )
    else:
        verdict["verdict_binding"] = True

    class_rows: list[dict[str, Any]] = []
    for model_id in sorted(model_scores):
        class_rows.extend(
            {"model_id": model_id, **row}
            for row in class_statistics(observed, model_scores[model_id])
        )
    seed_metric_rows: list[dict[str, Any]] = []
    for model_id in sorted(model_seed_scores):
        for index, seed in enumerate(model_seed_roster[model_id]):
            value = metrics(observed, model_seed_scores[model_id][index])
            seed_metric_rows.append(
                {
                    "model_id": model_id,
                    "model_seed": seed,
                    **{metric_id: value[metric_id] for metric_id in METRIC_IDS},
                }
            )

    write_tsv(output / "metrics.tsv", ("model_id", "metric", "estimate"), metric_rows)
    write_tsv(
        output / "bootstrap_intervals.tsv",
        ("model_id", "metric", "estimate", "ci_low", "ci_high", "bootstrap_replicates_used"),
        bootstrap_rows,
    )
    write_tsv(
        output / "paired_vs_comparator.tsv",
        tuple(paired_rows[0]),
        paired_rows,
    )
    write_tsv(output / "permutation_null.tsv", tuple(null_rows[0]), null_rows)
    write_tsv(
        output / "per_class.tsv",
        ("model_id", "class", "support", "precision", "recall", "f1", "tp", "fp", "fn"),
        class_rows,
    )
    write_tsv(
        output / "seed_level_metrics.tsv", tuple(seed_metric_rows[0]), seed_metric_rows
    )
    (output / "promotion_gate_verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "seed_direction_consistency.json").write_text(
        json.dumps(seed_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    receipt = {
        "schema_version": "masld-bench-gse49541-external-evaluation-v1",
        "status": "passed_external_development_scoring",
        "task_id": "gse49541_fibrosis_transfer",
        "arm_id": arm_id,
        "gate_eligible_arm": bool(gate_eligible_arm),
        "participants": EXPECTED_PARTICIPANTS,
        "biological_unit": "participant",
        "fibrosis_group_counts": dict(EXPECTED_GROUP_COUNTS),
        "external_prevalence_advanced_f3_f4": prevalence,
        "prevalence_is_not_the_auprc_null": True,
        "models": sorted(model_scores),
        "mandatory_baselines": list(MANDATORY_BASELINES),
        "unfittable_baselines_declared": declared,
        "unfittable_baseline_rationale": (
            "recorded by the source fit as having no registered configuration "
            "that converged and kept a non-zero coefficient across every outer "
            "fold and the full fit"
            if declared
            else "none"
        ),
        "scored_baselines": sorted(set(MANDATORY_BASELINES) & set(model_scores)),
        "candidate_model_id": candidate_model_id,
        "comparator_model_id": comparator_id,
        "primary_metric": "participant_auprc_advanced_f3_f4",
        "metrics_reported": list(METRIC_IDS),
        "brier_definition": "binary_brier_on_the_advanced_f3_f4_indicator",
        "bootstrap_method": "unstratified_participant_resample_paired_across_models",
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_replicates_used": int(np.sum(usable)),
        "bootstrap_seed": bootstrap_seed,
        "permutation_method": "outcome_permuted_against_each_models_own_realised_score_vector",
        "permutation_replicates": permutation_replicates,
        "permutation_seed": permutation_seed,
        "promotion_gate_verdict": verdict["verdict"],
        "promotion_gate_conditions_met": verdict["conditions_met"],
        "promotion_gate_conditions_total": verdict["conditions_total"],
        "query_fit_or_calibration_performed": False,
        "model_selected_from_external_outcomes": False,
        "external_outcomes_used_for_model_repair": False,
        "prediction_file_sha256": dict(sorted(prediction_file_sha256.items())),
        "prediction_bundle_artifacts_sha256": dict(
            sorted(prediction_artifacts_sha256.items())
        ),
        "project_sealed": False,
        "external_development_only": True,
        "champion_claim_eligible": False,
        "clinical_claim_eligible": False,
        "diagnostic_or_prognostic_claim_eligible": False,
        "confirmatory_claim_eligible": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluator-only", required=True, type=Path)
    parser.add_argument("--evaluator-only-artifacts-sha256", required=True)
    parser.add_argument("--gate", required=True, type=Path)
    parser.add_argument("--gate-sha256", required=True)
    parser.add_argument("--candidate-model-id", required=True)
    parser.add_argument("--arm-id", required=True)
    parser.add_argument("--gate-eligible-arm", required=True, choices=("true", "false"))
    parser.add_argument("--unfittable-baseline", action="append", default=[])
    parser.add_argument("--prediction", action="append", required=True)
    parser.add_argument("--prediction-artifacts-sha256", action="append", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    if (
        sha256_file(arguments.evaluator_only / "ARTIFACTS.json")
        != arguments.evaluator_only_artifacts_sha256
    ):
        raise ExternalTransferEvaluatorError("evaluator-only ARTIFACTS SHA differs")
    if sha256_file(arguments.gate) != arguments.gate_sha256:
        raise ExternalTransferEvaluatorError("promotion gate SHA differs")
    roots = [Path(value) for value in arguments.prediction]
    if len(roots) != len(arguments.prediction_artifacts_sha256):
        raise ExternalTransferEvaluatorError("prediction roots and SHA rosters differ")
    for root, expected in zip(roots, arguments.prediction_artifacts_sha256, strict=True):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise ExternalTransferEvaluatorError(f"prediction ARTIFACTS SHA differs: {root}")
    result = evaluate(
        labels_path=arguments.evaluator_only / "labels.tsv",
        prediction_roots=roots,
        gate_path=arguments.gate,
        candidate_model_id=arguments.candidate_model_id,
        output=arguments.output,
        arm_id=arguments.arm_id,
        gate_eligible_arm=arguments.gate_eligible_arm == "true",
        unfittable_baselines=arguments.unfittable_baseline,
        prediction_bundle_artifacts_sha256=arguments.prediction_artifacts_sha256,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
