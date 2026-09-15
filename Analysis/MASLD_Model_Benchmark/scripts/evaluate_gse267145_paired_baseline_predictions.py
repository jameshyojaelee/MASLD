#!/usr/bin/env python3
"""Independently evaluate held-participant GSE267145 H3 predictions/retrieval."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file


MODEL_IDS = (
    "training_mean_h3_profile",
    "pca_ridge",
    "reduced_rank_regression",
    "pls2",
    "sparse_cca",
)
PROFILE_MODEL_IDS = MODEL_IDS[:-1]
EXPECTED_PARTICIPANTS = 99
EXPECTED_H3_FEATURES = 96_460
EXPECTED_FOLD_COUNTS = {0: 21, 1: 21, 2: 21, 3: 19, 4: 17}
BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_SEED = 267145


class PairedBaselineEvaluationError(ValueError):
    """Raised when prediction or evaluator inputs violate the separation."""


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PairedBaselineEvaluationError(f"TSV lacks a header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
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


def multinomial_deviance(observed: np.ndarray, predicted: np.ndarray) -> float:
    counts = np.asarray(observed, dtype=np.float64)
    profile = np.asarray(predicted, dtype=np.float64)
    if (
        counts.ndim != 1
        or profile.shape != counts.shape
        or not np.all(np.isfinite(counts))
        or not np.all(np.isfinite(profile))
        or np.any(counts < 0.0)
        or np.any(profile <= 0.0)
        or not math.isclose(float(profile.sum()), 1.0, rel_tol=0.0, abs_tol=1.0e-6)
        or float(counts.sum()) <= 0.0
    ):
        raise PairedBaselineEvaluationError("deviance inputs differ")
    selected = counts > 0.0
    expected = float(counts.sum()) * profile[selected]
    value = float(
        2.0
        * np.sum(counts[selected] * (np.log(counts[selected]) - np.log(expected)))
    )
    if value < -1.0e-8:
        raise PairedBaselineEvaluationError("multinomial deviance is negative")
    return max(value, 0.0)


def tie_aware_retrieval(similarities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray(similarities, dtype=np.float64)
    if (
        matrix.ndim != 2
        or matrix.shape[0] != matrix.shape[1]
        or not np.all(np.isfinite(matrix))
    ):
        raise PairedBaselineEvaluationError("retrieval similarity matrix differs")
    reciprocal_rank = np.empty(len(matrix), dtype=np.float64)
    top1 = np.empty(len(matrix), dtype=np.float64)
    for index, row in enumerate(matrix):
        target = row[index]
        greater = int(np.sum(row > target + 1.0e-12))
        tied = int(np.sum(np.abs(row - target) <= 1.0e-12))
        if tied < 1:
            raise PairedBaselineEvaluationError("retrieval target is not in its tie set")
        average_rank = 1.0 + greater + 0.5 * (tied - 1)
        reciprocal_rank[index] = 1.0 / average_rank
        top1[index] = 1.0 / tied if greater == 0 else 0.0
    return reciprocal_rank, top1


def row_normalize(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or not np.all(np.isfinite(array)):
        raise PairedBaselineEvaluationError("retrieval embeddings differ")
    norms = np.linalg.norm(array, axis=1)
    if np.any(norms <= 0.0):
        raise PairedBaselineEvaluationError("retrieval embeddings differ")
    return array / norms[:, np.newaxis]


def fold_stratified_multiplicities(folds: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    result = np.zeros(
        (BOOTSTRAP_REPLICATES, EXPECTED_PARTICIPANTS), dtype=np.int16
    )
    for fold in range(5):
        indices = np.flatnonzero(folds == fold)
        if len(indices) != EXPECTED_FOLD_COUNTS[fold]:
            raise PairedBaselineEvaluationError("bootstrap fold census differs")
        draws = rng.integers(0, len(indices), size=(BOOTSTRAP_REPLICATES, len(indices)))
        row_axis = np.repeat(np.arange(BOOTSTRAP_REPLICATES), len(indices))
        column_axis = indices[draws.reshape(-1)]
        np.add.at(result, (row_axis, column_axis), 1)
    if np.any(result.sum(axis=1) != EXPECTED_PARTICIPANTS):
        raise PairedBaselineEvaluationError("bootstrap participant count differs")
    return result


def interval(values: np.ndarray) -> dict[str, Any]:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (BOOTSTRAP_REPLICATES,) or not np.all(np.isfinite(array)):
        raise PairedBaselineEvaluationError("bootstrap distribution differs")
    return {
        "lower": float(np.quantile(array, 0.025, method="linear")),
        "median": float(np.quantile(array, 0.5, method="linear")),
        "upper": float(np.quantile(array, 0.975, method="linear")),
        "confidence_level": 0.95,
        "replicates": BOOTSTRAP_REPLICATES,
        "seed": BOOTSTRAP_SEED,
    }


def run(
    *,
    evaluator_view: Path,
    predictions: Path,
    task_spec: Path,
    promotion_gate: Path,
    output: Path,
    expected_evaluator_artifacts_sha256: str,
    expected_prediction_artifacts_sha256: str,
    expected_task_spec_sha256: str,
    expected_promotion_gate_sha256: str,
) -> None:
    if output.exists():
        raise PairedBaselineEvaluationError("refusing to overwrite paired evaluation")
    for root, expected in (
        (evaluator_view, expected_evaluator_artifacts_sha256),
        (predictions, expected_prediction_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise PairedBaselineEvaluationError("input ARTIFACTS SHA-256 differs")
        verify_frozen_tree(root)
    if sha256_file(task_spec) != expected_task_spec_sha256:
        raise PairedBaselineEvaluationError("TaskSpec SHA-256 differs")
    if sha256_file(promotion_gate) != expected_promotion_gate_sha256:
        raise PairedBaselineEvaluationError("promotion gate SHA-256 differs")
    gate = json.loads(promotion_gate.read_text(encoding="utf-8"))
    if gate.get("champion_eligible") is not False or gate.get("claim_mode") != "development_only":
        raise PairedBaselineEvaluationError("promotion gate is not development-only")

    evaluator_metadata = verify_frozen_tree(evaluator_view)["metadata"]
    prediction_metadata = verify_frozen_tree(predictions)["metadata"]
    evaluator_receipt = json.loads(
        (evaluator_view / "receipt.json").read_text(encoding="utf-8")
    )
    prediction_receipt = json.loads(
        (predictions / "prediction_receipt.json").read_text(encoding="utf-8")
    )
    if (
        evaluator_metadata.get("artifact_class")
        != "gse267145_paired_h3_evaluator_view"
        or evaluator_metadata.get("participants") != EXPECTED_PARTICIPANTS
        or evaluator_metadata.get("h3k27ac_opaque_features") != EXPECTED_H3_FEATURES
        or evaluator_metadata.get("rna_values_included") is not False
        or evaluator_metadata.get("histology_or_participant_covariates_included")
        is not False
        or evaluator_receipt.get("h3_library_totals_evaluator_only") is not True
        or evaluator_receipt.get("coordinate_semantics_inferred") is not False
        or evaluator_receipt.get("task_spec_sha256") != expected_task_spec_sha256
        or evaluator_receipt.get("promotion_gate_sha256")
        != expected_promotion_gate_sha256
        or prediction_metadata.get("artifact_class")
        != "gse267145_paired_baseline_prediction_campaign"
        or prediction_metadata.get("outer_folds") != 5
        or prediction_metadata.get("held_participant_h3_read") is not False
        or prediction_metadata.get("metrics_calculated") is not False
        or prediction_receipt.get("model_ids") != list(MODEL_IDS)
        or prediction_receipt.get("folds_complete") is not True
        or prediction_receipt.get("task_spec_sha256") != expected_task_spec_sha256
        or prediction_receipt.get("promotion_gate_sha256")
        != expected_promotion_gate_sha256
        or prediction_receipt.get("champion_claim_allowed") is not False
    ):
        raise PairedBaselineEvaluationError("evaluator or prediction receipt differs")

    participant_fields, participant_rows = read_tsv(
        evaluator_view / "participant_axis.tsv"
    )
    if participant_fields != [
        "participant_index",
        "participant_id",
        "outer_fold",
        "h3k27ac_observation_state",
    ]:
        raise PairedBaselineEvaluationError("evaluator participant schema differs")
    participants = [row["participant_id"] for row in participant_rows]
    folds = np.asarray([int(row["outer_fold"]) for row in participant_rows], dtype=np.int8)
    if (
        len(participants) != EXPECTED_PARTICIPANTS
        or len(set(participants)) != EXPECTED_PARTICIPANTS
        or any(row["h3k27ac_observation_state"] != "observed" for row in participant_rows)
        or {fold: int(np.sum(folds == fold)) for fold in range(5)}
        != EXPECTED_FOLD_COUNTS
    ):
        raise PairedBaselineEvaluationError("evaluator participant axis differs")
    observed = np.load(
        evaluator_view / "observed_h3k27ac_counts.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    if observed.shape != (EXPECTED_PARTICIPANTS, EXPECTED_H3_FEATURES) or observed.dtype != np.uint32:
        raise PairedBaselineEvaluationError("observed evaluator H3 array differs")
    observed_totals = observed.sum(axis=1, dtype=np.float64)
    if np.any(observed_totals <= 0.0):
        raise PairedBaselineEvaluationError("observed evaluator H3 total is nonpositive")

    profile_deviance = {
        model_id: np.full(EXPECTED_PARTICIPANTS, np.nan, dtype=np.float64)
        for model_id in PROFILE_MODEL_IDS
    }
    retrieval_rr = {
        model_id: np.full(EXPECTED_PARTICIPANTS, np.nan, dtype=np.float64)
        for model_id in MODEL_IDS
    }
    retrieval_top1 = {
        model_id: np.full(EXPECTED_PARTICIPANTS, np.nan, dtype=np.float64)
        for model_id in MODEL_IDS
    }
    participant_to_index = {participant: index for index, participant in enumerate(participants)}

    for outer_fold in range(5):
        fold_root = predictions / f"outer_{outer_fold}"
        verify_frozen_tree(fold_root)
        fold_receipt = json.loads(
            (fold_root / "prediction_receipt.json").read_text(encoding="utf-8")
        )
        query_fields, query_rows = read_tsv(fold_root / "query_participants.tsv")
        if query_fields != ["participant_index", "participant_id", "outer_fold"]:
            raise PairedBaselineEvaluationError("prediction participant schema differs")
        query_indices = np.asarray(
            [participant_to_index[row["participant_id"]] for row in query_rows],
            dtype=np.int16,
        )
        if (
            len(query_rows) != EXPECTED_FOLD_COUNTS[outer_fold]
            or query_indices.tolist()
            != [int(row["participant_index"]) for row in query_rows]
            or np.any(folds[query_indices] != outer_fold)
            or fold_receipt.get("outer_fold") != outer_fold
            or fold_receipt.get("model_ids") != list(MODEL_IDS)
            or fold_receipt.get("held_participant_h3_read") is not False
            or fold_receipt.get("histology_sex_age_nas_fibrosis_read") is not False
            or fold_receipt.get("metrics_calculated") is not False
            or fold_receipt.get("sparse_cca_profile_prediction_applicable") is not False
            or fold_receipt.get("sparse_cca_retrieval_applicable") is not True
        ):
            raise PairedBaselineEvaluationError("fold prediction receipt differs")
        fold_counts = np.asarray(observed[query_indices], dtype=np.uint32)
        fold_totals = observed_totals[query_indices]
        fold_profiles = fold_counts.astype(np.float64) / fold_totals[:, np.newaxis]
        observed_hellinger = np.sqrt(fold_profiles)

        for model_id in PROFILE_MODEL_IDS:
            predicted = np.load(
                fold_root / "predictions" / f"{model_id}.npy",
                allow_pickle=False,
            )
            if (
                predicted.shape != (len(query_indices), EXPECTED_H3_FEATURES)
                or predicted.dtype != np.float32
                or not np.all(np.isfinite(predicted))
                or np.any(predicted <= 0.0)
                or not np.allclose(
                    predicted.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-6
                )
            ):
                raise PairedBaselineEvaluationError("profile prediction array differs")
            for local, global_index in enumerate(query_indices):
                profile_deviance[model_id][global_index] = multinomial_deviance(
                    fold_counts[local], predicted[local]
                )
            similarities = np.sqrt(predicted.astype(np.float64)) @ observed_hellinger.T
            reciprocal, top1 = tie_aware_retrieval(similarities)
            retrieval_rr[model_id][query_indices] = reciprocal
            retrieval_top1[model_id][query_indices] = top1

        query_embedding = np.load(
            fold_root / "sparse_cca/query_embeddings.npy", allow_pickle=False
        )
        with np.load(
            fold_root / "sparse_cca/candidate_transform.npz", allow_pickle=False
        ) as state:
            required = {
                "h3_feature_indices",
                "h3_mean",
                "h3_scale",
                "h3_weights",
            }
            if set(state.files) != required:
                raise PairedBaselineEvaluationError("sparse CCA candidate state differs")
            feature_indices = state["h3_feature_indices"]
            h3_mean = state["h3_mean"]
            h3_scale = state["h3_scale"]
            h3_weights = state["h3_weights"]
        if (
            query_embedding.shape != (len(query_indices), h3_weights.shape[1])
            or feature_indices.ndim != 1
            or len(feature_indices) != len(h3_mean)
            or len(feature_indices) != len(h3_scale)
            or h3_weights.shape[0] != len(feature_indices)
            or len(set(feature_indices.tolist())) != len(feature_indices)
            or np.any(feature_indices < 0)
            or np.any(feature_indices >= EXPECTED_H3_FEATURES)
            or np.any(h3_scale <= 0.0)
        ):
            raise PairedBaselineEvaluationError("sparse CCA transform dimensions differ")
        candidate = (
            observed_hellinger[:, feature_indices] - h3_mean.astype(np.float64)
        ) / h3_scale.astype(np.float64)
        candidate_embedding = candidate @ h3_weights.astype(np.float64)
        similarities = row_normalize(query_embedding) @ row_normalize(candidate_embedding).T
        reciprocal, top1 = tie_aware_retrieval(similarities)
        retrieval_rr["sparse_cca"][query_indices] = reciprocal
        retrieval_top1["sparse_cca"][query_indices] = top1

    if any(not np.all(np.isfinite(values)) for values in profile_deviance.values()):
        raise PairedBaselineEvaluationError("profile evaluation coverage differs")
    if any(not np.all(np.isfinite(values)) for values in retrieval_rr.values()):
        raise PairedBaselineEvaluationError("retrieval evaluation coverage differs")

    multiplicities = fold_stratified_multiplicities(folds).astype(np.float64)
    denominators = multiplicities.sum(axis=1)
    baseline_deviance = profile_deviance["training_mean_h3_profile"]
    if np.any(baseline_deviance <= 0.0):
        raise PairedBaselineEvaluationError(
            "training-mean participant deviance is nonpositive"
        )
    profile_rows: list[dict[str, Any]] = []
    participant_profile_rows: list[dict[str, Any]] = []
    retrieval_rows: list[dict[str, Any]] = []
    participant_retrieval_rows: list[dict[str, Any]] = []
    interval_rows: list[dict[str, Any]] = []
    distributions: dict[str, np.ndarray] = {}

    for model_id in PROFILE_MODEL_IDS:
        values = profile_deviance[model_id]
        mean_deviance = float(np.mean(values))
        participant_skill = 1.0 - values / baseline_deviance
        skill = float(np.mean(participant_skill))
        aggregate_skill = 1.0 - mean_deviance / float(np.mean(baseline_deviance))
        profile_rows.append(
            {
                "model_id": model_id,
                "participants": EXPECTED_PARTICIPANTS,
                "participant_macro_multinomial_deviance": mean_deviance,
                "participant_macro_deviance_skill_vs_training_mean": skill,
                "aggregate_deviance_skill_depth_sensitive_diagnostic": aggregate_skill,
            }
        )
        for index, participant in enumerate(participants):
            participant_profile_rows.append(
                {
                    "model_id": model_id,
                    "participant_id": participant,
                    "outer_fold": int(folds[index]),
                    "observed_h3_library_total": int(observed_totals[index]),
                    "multinomial_deviance": values[index],
                    "training_mean_multinomial_deviance": baseline_deviance[index],
                    "deviance_skill_vs_training_mean": participant_skill[index],
                }
            )
        boot_deviance = multiplicities @ values / denominators
        boot_skill = multiplicities @ participant_skill / denominators
        for metric, observed_value, distribution in (
            ("participant_macro_multinomial_deviance", mean_deviance, boot_deviance),
            (
                "participant_macro_deviance_skill_vs_training_mean",
                skill,
                boot_skill,
            ),
        ):
            distributions[f"{model_id}::{metric}"] = distribution
            interval_rows.append(
                {
                    "lane": "rna_conditioned_h3_profile",
                    "model_id": model_id,
                    "metric": metric,
                    "estimate": observed_value,
                    **interval(distribution),
                }
            )

    for model_id in MODEL_IDS:
        rr = retrieval_rr[model_id]
        top1 = retrieval_top1[model_id]
        mrr = float(np.mean(rr))
        top1_mean = float(np.mean(top1))
        retrieval_rows.append(
            {
                "model_id": model_id,
                "participants": EXPECTED_PARTICIPANTS,
                "paired_participant_retrieval_mrr": mrr,
                "paired_participant_retrieval_top1": top1_mean,
            }
        )
        for index, participant in enumerate(participants):
            participant_retrieval_rows.append(
                {
                    "model_id": model_id,
                    "participant_id": participant,
                    "outer_fold": int(folds[index]),
                    "candidate_pool_size": EXPECTED_FOLD_COUNTS[int(folds[index])],
                    "reciprocal_rank": rr[index],
                    "top1_tie_credit": top1[index],
                }
            )
        for metric, observed_value, values in (
            ("paired_participant_retrieval_mrr", mrr, rr),
            ("paired_participant_retrieval_top1", top1_mean, top1),
        ):
            distribution = multiplicities @ values / denominators
            distributions[f"{model_id}::{metric}"] = distribution
            interval_rows.append(
                {
                    "lane": "observed_pair_retrieval",
                    "model_id": model_id,
                    "metric": metric,
                    "estimate": observed_value,
                    **interval(distribution),
                }
            )

    output.mkdir(mode=0o750)
    write_tsv(output / "model_profile_metrics.tsv", list(profile_rows[0]), profile_rows)
    write_tsv(
        output / "participant_profile_metrics.tsv",
        list(participant_profile_rows[0]),
        participant_profile_rows,
    )
    write_tsv(
        output / "model_retrieval_metrics.tsv", list(retrieval_rows[0]), retrieval_rows
    )
    write_tsv(
        output / "participant_retrieval_metrics.tsv",
        list(participant_retrieval_rows[0]),
        participant_retrieval_rows,
    )
    write_tsv(output / "bootstrap_intervals.tsv", list(interval_rows[0]), interval_rows)
    dispositions = [
        {
            "model_id": model_id,
            "rna_conditioned_h3_profile": (
                "evaluated" if model_id in PROFILE_MODEL_IDS else "not_applicable"
            ),
            "observed_pair_retrieval": "evaluated",
            "reason": (
                "directional_profile_prediction_available"
                if model_id in PROFILE_MODEL_IDS
                else "symmetric_sparse_cca_not_forced_into_directional_profile_regression"
            ),
        }
        for model_id in MODEL_IDS
    ]
    write_tsv(output / "model_endpoint_dispositions.tsv", list(dispositions[0]), dispositions)
    with (output / "bootstrap_distributions.npz").open("xb") as handle:
        np.savez_compressed(handle, **distributions)
    receipt = {
        "schema_version": "masld-bench-gse267145-paired-independent-evaluator-v1",
        "status": "pass_development_evaluation_only",
        "task_id": "paired_bulk_rna_h3k27ac",
        "participants": EXPECTED_PARTICIPANTS,
        "biological_unit": "participant",
        "pairing": "same_sample_different_aliquot",
        "model_ids": list(MODEL_IDS),
        "profile_models": list(PROFILE_MODEL_IDS),
        "retrieval_models": list(MODEL_IDS),
        "profile_endpoint": {
            "metric": "participant_macro_h3_profile_multinomial_deviance_skill",
            "reference": "outer_training_mean_h3_profile",
            "held_participant_library_total_read_only_inside_evaluator": True,
            "participants_receive_equal_mass": True,
            "skill_definition": "arithmetic_mean_across_participants_of_one_minus_model_deviance_over_that_participants_training_mean_deviance",
        },
        "retrieval_endpoint": {
            "metrics": [
                "paired_participant_retrieval_mrr",
                "paired_participant_retrieval_top1",
            ],
            "candidate_pool": "held_participants_within_same_outer_fold",
            "profile_similarity": "hellinger_bhattacharyya_coefficient",
            "sparse_cca_similarity": "cosine_in_training_fitted_canonical_space",
            "ties": "average_rank_and_fractional_top1_credit",
            "bootstrap_conditions_on_frozen_candidate_pools": True,
            "cannot_support_rna_conditioned_profile_claim": True,
        },
        "bootstrap": {
            "method": "outer_fold_stratified_paired_participant_cluster_bootstrap",
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": BOOTSTRAP_SEED,
            "participants_are_independent_units": True,
            "profiles_or_features_are_independent_units": False,
            "distribution_sha256": canonical_sha256(
                {key: value.tolist() for key, value in sorted(distributions.items())}
            ),
        },
        "task_spec_sha256": expected_task_spec_sha256,
        "promotion_gate_sha256": expected_promotion_gate_sha256,
        "evaluator_view_artifacts_sha256": expected_evaluator_artifacts_sha256,
        "prediction_artifacts_sha256": expected_prediction_artifacts_sha256,
        "models_fit_inside_evaluator": False,
        "preprocessing_fit_inside_evaluator": False,
        "held_participant_h3_read_only_inside_evaluator": True,
        "histology_sex_age_nas_fibrosis_read": False,
        "coordinate_interpretation_used": False,
        "multiplicity_family": "paired_bulk_rna_h3k27ac_development_secondary",
        "multiplicity_tests_performed": False,
        "model_selection_performed": False,
        "champion_promotion_performed": False,
        "external_evaluation_performed": False,
        "artifact_release_state": "internal_only_source_data_redistribution_prohibited",
        "clinical_claim_allowed": False,
    }
    write_json_exclusive(output / "evaluation_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse267145_paired_baseline_development_evaluation",
            "participants": EXPECTED_PARTICIPANTS,
            "models": len(MODEL_IDS),
            "models_fit_inside_evaluator": False,
            "model_selection_performed": False,
            "champion_promotion_performed": False,
            "external_evaluation_performed": False,
            "status": "passed",
        },
    )


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--evaluator-view", required=True, type=Path)
    value.add_argument("--evaluator-artifacts-sha256", required=True)
    value.add_argument("--predictions", required=True, type=Path)
    value.add_argument("--prediction-artifacts-sha256", required=True)
    value.add_argument("--task-spec", required=True, type=Path)
    value.add_argument("--task-spec-sha256", required=True)
    value.add_argument("--promotion-gate", required=True, type=Path)
    value.add_argument("--promotion-gate-sha256", required=True)
    value.add_argument("--output", required=True, type=Path)
    return value


def main() -> int:
    arguments = parser().parse_args()
    run(
        evaluator_view=arguments.evaluator_view,
        predictions=arguments.predictions,
        task_spec=arguments.task_spec,
        promotion_gate=arguments.promotion_gate,
        output=arguments.output,
        expected_evaluator_artifacts_sha256=arguments.evaluator_artifacts_sha256,
        expected_prediction_artifacts_sha256=arguments.prediction_artifacts_sha256,
        expected_task_spec_sha256=arguments.task_spec_sha256,
        expected_promotion_gate_sha256=arguments.promotion_gate_sha256,
    )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
