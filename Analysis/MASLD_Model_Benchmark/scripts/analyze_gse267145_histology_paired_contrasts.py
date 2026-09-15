#!/usr/bin/env python3
"""Paired development contrasts over the frozen GSE267145 v4 histology scores.

No model is refitted and no frozen prediction bundle is mutated.  The frozen
scoring campaign's participant bootstrap is regenerated deterministically from
its recorded seed, verified against its recorded digest, and reused so that
every contrast is paired on the same resamples.

The metric engine is the frozen scorer's own ``_bootstrap_metrics``: the whole
metric is recomputed on each resampled participant set for both arms and then
differenced.  That is required because five of the nine metrics are set-level
statistics with no per-participant decomposition.

No p-value and no BH adjustment is produced.  The requirements at
``config/evaluation/gse267145_histology_production_v4_scoring.json`` sets
``confirmatory_inference_allowed`` to false, and the adoption criteria record the
multiplicity family as not applicable without a confirmatory claim.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file

import scripts.evaluate_gse267145_histology_predictions as reference
import scripts.score_gse267145_histology_production_v4 as scorer


PRESPEC_SHA256_FIELD = "prespec_contract_sha256"
ANALYSIS_REVISION_ID = "model-check-1102-paired-contrast-v1"
EXPECTED_PARTICIPANTS = 99
BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_SEED = 20260824
BOOTSTRAP_INDICES_SHA256 = (
    "fb7ea4b28ef37ee759a70bb3664868f6affa5accc8dadd96edf202238536b2c5"
)
SCORER_SOURCE_SHA256 = (
    "dad1fa5b68686eaea6eeaa9c9725a4b0aed23dd4e96fcedbe1857d02ea070a43"
)
REFERENCE_EVALUATOR_SHA256 = (
    "1d753b8e23fda8f9582e4d8cc2ce9c90eeae65d370ce342bf9d32ffdd2ba5920"
)
SCORES_ARTIFACTS_SHA256 = (
    "a55aed3d9ac4cad2211ed3b319c45f911d05851b8559408ed99bf7fa59562393"
)
OUTCOMES_ARTIFACTS_SHA256 = (
    "98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a"
)
FOLDS_ARTIFACTS_SHA256 = (
    "9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3"
)

METRIC_DIRECTION = {
    "stage3_macro_f1": 1.0,
    "stage3_multiclass_brier": -1.0,
    "fibrosis_group3_macro_f1": 1.0,
    "fibrosis_cumulative_ordinal_mae": -1.0,
    "fibrosis_cumulative_spearman": 1.0,
    "fibrosis_regression_mae": -1.0,
    "fibrosis_regression_spearman": 1.0,
    "nash_crn_component_sum_spearman": 1.0,
    "nash_crn_component_sum_mae": -1.0,
}
METRIC_IDS = tuple(METRIC_DIRECTION)

MATCHED_MODALITY_PAIRS = (
    ("h3_variance_pca_elastic_net", "rna_hvg_pca_elastic_net", "elastic_net"),
    ("h3_variance_pca_linear_svm", "rna_hvg_pca_linear_svm", "linear_svm"),
    ("h3_variance_pca_nearest_centroid", "rna_hvg_pca_nearest_centroid", "nearest_centroid"),
    ("h3_variance_pca_knn", "rna_hvg_pca_knn", "knn"),
)
FUSION_ARMS = ("block_pca_elastic_net", "calibrated_late_fusion")
FUSION_PARENTS = ("rna_hvg_pca_elastic_net", "h3_variance_pca_elastic_net")
BASELINE_MODEL_ID = "training_stage_distribution"
BASELINE_CONTRAST_ROSTER = (
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

SUMMARY_FIELDS = (
    "contrast_id",
    "family_id",
    "metric_id",
    "direction",
    "arm_a",
    "arm_b",
    "applicability_state",
    "applicability_reason",
    "arm_a_estimate",
    "arm_b_estimate",
    "signed_difference",
    "ci95_low",
    "ci95_high",
    "probability_improvement",
    "valid_bootstrap_replicates",
    "interval_excludes_zero",
    "biological_resampling_unit",
    "p_value",
    "bh_adjusted_q_value",
    "confirmatory_inference_allowed",
    "claim_mode",
)


class PairedContrastError(RuntimeError):
    """Raised when a bound input or a held-back criterion does not hold."""


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _verify_file(path: Path, expected: str, label: str) -> None:
    if not path.is_file() or path.is_symlink() or sha256_file(path) != expected:
        raise PairedContrastError(f"{label} SHA-256 differs")


def _format(value: float) -> str:
    return "not_estimable" if not math.isfinite(value) else format(value, ".17g")


def _bootstrap_indices() -> np.ndarray:
    indices = np.random.default_rng(BOOTSTRAP_SEED).integers(
        0,
        EXPECTED_PARTICIPANTS,
        size=(BOOTSTRAP_REPLICATES, EXPECTED_PARTICIPANTS),
        endpoint=False,
    )
    if canonical_sha256(indices.tolist()) != BOOTSTRAP_INDICES_SHA256:
        raise PairedContrastError(
            "regenerated bootstrap indices do not match the frozen scoring digest"
        )
    return indices


def _participant_axis(folds: Path) -> tuple[list[str], list[str]]:
    fields, rows = scorer._read_tsv(folds / "participant_outer_folds.tsv")
    if fields != ("participant_id", "outer_fold") or len(rows) != EXPECTED_PARTICIPANTS:
        raise PairedContrastError("participant fold axis differs")
    return (
        [row["participant_id"] for row in rows],
        [row["outer_fold"] for row in rows],
    )


def _load_model_predictions(
    bundle_root: Path, participants: Sequence[str], outer_folds: Sequence[str]
) -> dict[str, dict[str, np.ndarray]]:
    """Reassemble each model's five frozen endpoint bundles into one array set."""

    verify_frozen_tree(bundle_root)
    result: dict[str, dict[str, np.ndarray]] = {}
    for model_id in reference.MODEL_IDS:
        values: dict[str, np.ndarray] = {}
        for endpoint_id in scorer.PREDICTED_ENDPOINTS:
            target = bundle_root / f"{model_id}--endpoint-{endpoint_id}"
            verify_frozen_tree(target)
            _fields, rows = scorer._read_tsv(target / "predictions.tsv")
            if (
                len(rows) != EXPECTED_PARTICIPANTS
                or [row["participant_id"] for row in rows] != list(participants)
                or [row["outer_fold"] for row in rows] != list(outer_folds)
                or any(row["endpoint_id"] != endpoint_id for row in rows)
                or any(row["prediction_state"] != "observed" for row in rows)
            ):
                raise PairedContrastError(
                    f"frozen prediction bundle axis differs: {model_id}/{endpoint_id}"
                )
            if endpoint_id == "stage3":
                values["stage_probabilities"] = np.asarray(
                    [
                        [float(row[f"probability_{label}"]) for label in reference.STAGE3]
                        for row in rows
                    ],
                    dtype=np.float64,
                )
                values["stage3"] = np.asarray(
                    [row["predicted_stage3"] for row in rows], dtype=str
                )
            elif endpoint_id == "fibrosis_group3":
                values["fibrosis_probabilities"] = np.asarray(
                    [
                        [
                            float(row[f"probability_fibrosis_{label}"])
                            for label in reference.FIBROSIS_GROUP3
                        ]
                        for row in rows
                    ],
                    dtype=np.float64,
                )
                values["fibrosis_group3"] = np.asarray(
                    [row["predicted_fibrosis_group3"] for row in rows], dtype=str
                )
            elif endpoint_id == "fibrosis_cumulative":
                values["fibrosis_cumulative"] = np.asarray(
                    [float(row["predicted_fibrosis_cumulative_expected"]) for row in rows],
                    dtype=np.float64,
                )
            elif endpoint_id == "fibrosis_regression":
                values["fibrosis_regression"] = np.asarray(
                    [float(row["predicted_fibrosis_regression"]) for row in rows],
                    dtype=np.float64,
                )
            else:
                values["nas"] = np.asarray(
                    [float(row["predicted_nash_crn_component_sum"]) for row in rows],
                    dtype=np.float64,
                )
        result[model_id] = scorer._validate_prediction_arrays(values)
    if len(result) != 11:
        raise PairedContrastError("frozen prediction bundle model roster differs")
    return result


def _frozen_point_metrics(scores: Path) -> dict[tuple[str, str], float]:
    """Read the frozen point estimates so the recomputation can be checked."""

    _fields, rows = scorer._read_tsv(scores / "standardized_metrics.tsv")
    frozen: dict[tuple[str, str], float] = {}
    for row in rows:
        if row["applicability_state"] != "observed":
            continue
        estimate = row["estimate"]
        frozen[(row["model_id"], row["metric_id"])] = (
            math.nan if estimate == "not_estimable" else float(estimate)
        )
    if len(frozen) != 99:
        raise PairedContrastError("frozen endpoint-native metric row count differs")
    return frozen


def _check_recomputation(
    points: Mapping[str, Mapping[str, float]],
    frozen: Mapping[tuple[str, str], float],
) -> int:
    checked = 0
    for (model_id, metric_id), frozen_value in frozen.items():
        recomputed = points[model_id][metric_id]
        if math.isnan(frozen_value):
            if not math.isnan(recomputed):
                raise PairedContrastError(
                    f"recomputed {model_id}/{metric_id} is finite but the frozen value is not"
                )
        elif not math.isclose(recomputed, frozen_value, rel_tol=0.0, abs_tol=1e-12):
            raise PairedContrastError(
                f"recomputed {model_id}/{metric_id} differs from the frozen point estimate"
            )
        checked += 1
    return checked


def _contrast(
    *,
    contrast_id: str,
    family_id: str,
    metric_id: str,
    arm_a: str,
    arm_b: str,
    arm_a_point: float,
    arm_b_point: float,
    arm_a_draws: np.ndarray,
    arm_b_draws: np.ndarray,
) -> tuple[dict[str, Any], np.ndarray | None]:
    direction = METRIC_DIRECTION[metric_id]
    row: dict[str, Any] = {
        "contrast_id": contrast_id,
        "family_id": family_id,
        "metric_id": metric_id,
        "direction": "maximize" if direction > 0 else "minimize",
        "arm_a": arm_a,
        "arm_b": arm_b,
        "arm_a_estimate": _format(arm_a_point),
        "arm_b_estimate": _format(arm_b_point),
        "biological_resampling_unit": "participant",
        "p_value": "not_calculated_development_no_confirmatory_lock",
        "bh_adjusted_q_value": "not_calculated_development_no_confirmatory_lock",
        "confirmatory_inference_allowed": "false",
        "claim_mode": "single_cohort_participant_held_development_only",
    }
    if not math.isfinite(arm_a_point) or not math.isfinite(arm_b_point):
        absent = arm_b if not math.isfinite(arm_b_point) else arm_a
        row.update(
            {
                "applicability_state": "not_applicable",
                "applicability_reason": (
                    "not_applicable_baseline_not_estimable"
                    if absent == BASELINE_MODEL_ID
                    else "not_applicable_arm_not_estimable"
                ),
                "signed_difference": "not_applicable",
                "ci95_low": "not_applicable",
                "ci95_high": "not_applicable",
                "probability_improvement": "not_applicable",
                "valid_bootstrap_replicates": 0,
                "interval_excludes_zero": "not_applicable",
            }
        )
        return row, None

    differences = direction * (arm_a_draws - arm_b_draws)
    valid = np.isfinite(differences)
    valid_count = int(np.sum(valid))
    if valid_count == 0:
        row.update(
            {
                "applicability_state": "not_applicable",
                "applicability_reason": "not_applicable_no_valid_replicate",
                "signed_difference": _format(direction * (arm_a_point - arm_b_point)),
                "ci95_low": "not_applicable",
                "ci95_high": "not_applicable",
                "probability_improvement": "not_applicable",
                "valid_bootstrap_replicates": 0,
                "interval_excludes_zero": "not_applicable",
            }
        )
        return row, None

    finite = differences[valid]
    low = float(np.percentile(finite, 2.5))
    high = float(np.percentile(finite, 97.5))
    row.update(
        {
            "applicability_state": "observed",
            "applicability_reason": "endpoint_native_paired_contrast",
            "signed_difference": _format(direction * (arm_a_point - arm_b_point)),
            "ci95_low": _format(low),
            "ci95_high": _format(high),
            "probability_improvement": _format(float(np.mean(finite > 0.0))),
            "valid_bootstrap_replicates": valid_count,
            "interval_excludes_zero": "true" if low > 0.0 or high < 0.0 else "false",
        }
    )
    return row, differences


def _build_contrasts(
    points: Mapping[str, Mapping[str, float]],
    draws: Mapping[str, Mapping[str, np.ndarray]],
) -> tuple[list[dict[str, Any]], dict[str, np.ndarray]]:
    rows: list[dict[str, Any]] = []
    vectors: dict[str, np.ndarray] = {}

    def emit(**kwargs: Any) -> None:
        row, difference = _contrast(**kwargs)
        rows.append(row)
        if difference is not None:
            vectors[row["contrast_id"]] = difference.astype(np.float64)

    for metric_id in METRIC_IDS:
        for h3_id, rna_id, algorithm in MATCHED_MODALITY_PAIRS:
            emit(
                contrast_id=f"M--{algorithm}--{metric_id}",
                family_id="M",
                metric_id=metric_id,
                arm_a=h3_id,
                arm_b=rna_id,
                arm_a_point=points[h3_id][metric_id],
                arm_b_point=points[rna_id][metric_id],
                arm_a_draws=draws[h3_id][metric_id],
                arm_b_draws=draws[rna_id][metric_id],
            )

    direction_by_metric = METRIC_DIRECTION
    for metric_id in METRIC_IDS:
        direction = direction_by_metric[metric_id]
        stacked = np.stack(
            [
                direction * (draws[h3_id][metric_id] - draws[rna_id][metric_id])
                for h3_id, rna_id, _ in MATCHED_MODALITY_PAIRS
            ]
        )
        mean_difference = np.mean(stacked, axis=0)
        point = float(
            np.mean(
                [
                    direction * (points[h3_id][metric_id] - points[rna_id][metric_id])
                    for h3_id, rna_id, _ in MATCHED_MODALITY_PAIRS
                ]
            )
        )
        contrast_id = f"MM--mean_of_four_matched_pairs--{metric_id}"
        valid = np.isfinite(mean_difference)
        valid_count = int(np.sum(valid))
        finite = mean_difference[valid]
        low = float(np.percentile(finite, 2.5)) if valid_count else math.nan
        high = float(np.percentile(finite, 97.5)) if valid_count else math.nan
        rows.append(
            {
                "contrast_id": contrast_id,
                "family_id": "MM",
                "metric_id": metric_id,
                "direction": "maximize" if direction > 0 else "minimize",
                "arm_a": "h3k27ac_only_lane_mean_of_four_matched_pairs",
                "arm_b": "rna_only_lane_mean_of_four_matched_pairs",
                "applicability_state": "observed" if valid_count else "not_applicable",
                "applicability_reason": (
                    "prespecified_mean_across_four_matched_algorithm_pairs"
                    if valid_count
                    else "not_applicable_no_valid_replicate"
                ),
                "arm_a_estimate": "not_applicable",
                "arm_b_estimate": "not_applicable",
                "signed_difference": _format(point),
                "ci95_low": _format(low),
                "ci95_high": _format(high),
                "probability_improvement": (
                    _format(float(np.mean(finite > 0.0))) if valid_count else "not_applicable"
                ),
                "valid_bootstrap_replicates": valid_count,
                "interval_excludes_zero": (
                    ("true" if low > 0.0 or high < 0.0 else "false")
                    if valid_count
                    else "not_applicable"
                ),
                "biological_resampling_unit": "participant",
                "p_value": "not_calculated_development_no_confirmatory_lock",
                "bh_adjusted_q_value": "not_calculated_development_no_confirmatory_lock",
                "confirmatory_inference_allowed": "false",
                "claim_mode": "single_cohort_participant_held_development_only",
            }
        )
        if valid_count:
            vectors[contrast_id] = mean_difference.astype(np.float64)

    for metric_id in METRIC_IDS:
        direction = METRIC_DIRECTION[metric_id]
        for arm in FUSION_ARMS:
            rna_parent, h3_parent = FUSION_PARENTS
            best_parent = max(
                FUSION_PARENTS,
                key=lambda parent: direction * points[parent][metric_id],
            )
            for label, parent in (
                ("rna_parent", rna_parent),
                ("h3_parent", h3_parent),
                ("best_parent", best_parent),
            ):
                emit(
                    contrast_id=f"F--{arm}--{label}--{metric_id}",
                    family_id="F",
                    metric_id=metric_id,
                    arm_a=arm,
                    arm_b=parent,
                    arm_a_point=points[arm][metric_id],
                    arm_b_point=points[parent][metric_id],
                    arm_a_draws=draws[arm][metric_id],
                    arm_b_draws=draws[parent][metric_id],
                )

    for metric_id in METRIC_IDS:
        for model_id in BASELINE_CONTRAST_ROSTER:
            emit(
                contrast_id=f"B--{model_id}--{metric_id}",
                family_id="B",
                metric_id=metric_id,
                arm_a=model_id,
                arm_b=BASELINE_MODEL_ID,
                arm_a_point=points[model_id][metric_id],
                arm_b_point=points[BASELINE_MODEL_ID][metric_id],
                arm_a_draws=draws[model_id][metric_id],
                arm_b_draws=draws[BASELINE_MODEL_ID][metric_id],
            )

    return rows, vectors


def _write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    path.chmod(0o440)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prespec", type=Path, required=True)
    parser.add_argument("--prespec-sha256", required=True)
    parser.add_argument("--sealed-prespec", type=Path, required=True)
    parser.add_argument("--sealed-prespec-artifacts-sha256", required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--folds", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    _verify_file(arguments.prespec, arguments.prespec_sha256, "prespec contract")
    _verify_file(
        arguments.sealed_prespec / "ARTIFACTS.json",
        arguments.sealed_prespec_artifacts_sha256,
        "sealed prespec campaign",
    )
    verify_frozen_tree(arguments.sealed_prespec)
    prespec = _read_json(arguments.prespec)
    if (
        prespec.get("status") != "sealed_before_any_contrast_is_computed"
        or prespec["resampling"]["bootstrap_indices_sha256"] != BOOTSTRAP_INDICES_SHA256
        or prespec["estimator"]["p_values_calculated"] is not False
        or prespec["estimator"]["bh_adjustment_calculated"] is not False
        or prespec["claim_boundary"]["model_ranking_allowed"] is not False
        or prespec["claim_boundary"]["confirmatory_inference_allowed"] is not False
        or prespec["declared_contrast_totals"]["declared"] != 189
        or prespec["declared_contrast_totals"]["computable"] != 169
    ):
        raise PairedContrastError("sealed prespec differs")

    root = Path(__file__).resolve().parent.parent
    _verify_file(
        root / "scripts/score_gse267145_histology_production_v4.py",
        SCORER_SOURCE_SHA256,
        "frozen scorer source",
    )
    _verify_file(
        root / "scripts/evaluate_gse267145_histology_predictions.py",
        REFERENCE_EVALUATOR_SHA256,
        "frozen reference evaluator",
    )
    _verify_file(
        arguments.scores / "ARTIFACTS.json", SCORES_ARTIFACTS_SHA256, "frozen scores"
    )
    _verify_file(
        arguments.outcomes / "ARTIFACTS.json", OUTCOMES_ARTIFACTS_SHA256, "outcomes"
    )
    _verify_file(arguments.folds / "ARTIFACTS.json", FOLDS_ARTIFACTS_SHA256, "folds")
    verify_frozen_tree(arguments.scores)
    verify_frozen_tree(arguments.outcomes)
    verify_frozen_tree(arguments.folds)

    indices = _bootstrap_indices()
    participants, outer_folds = _participant_axis(arguments.folds)
    outcomes = scorer._load_outcomes(arguments.outcomes, participants, outer_folds)
    predictions = _load_model_predictions(
        arguments.scores / "prediction_bundles", participants, outer_folds
    )

    points: dict[str, dict[str, float]] = {}
    draws: dict[str, dict[str, np.ndarray]] = {}
    for model_id in reference.MODEL_IDS:
        points[model_id] = scorer._point_metrics(outcomes, predictions[model_id])
        draws[model_id] = scorer._bootstrap_metrics(
            outcomes, predictions[model_id], indices
        )

    frozen = _frozen_point_metrics(arguments.scores)
    reproduced = _check_recomputation(points, frozen)

    rows, vectors = _build_contrasts(points, draws)
    observed = [row for row in rows if row["applicability_state"] == "observed"]
    not_applicable = [row for row in rows if row["applicability_state"] != "observed"]
    if len(rows) != 189 or len(observed) != 169 or len(not_applicable) != 20:
        raise PairedContrastError(
            f"contrast census differs: {len(rows)} declared, {len(observed)} observed"
        )

    arguments.output.mkdir(parents=True)
    _write_tsv(arguments.output / "paired_contrasts.tsv", SUMMARY_FIELDS, rows)
    np.savez_compressed(
        arguments.output / "paired_difference_draws.npz",
        **{key: value for key, value in sorted(vectors.items())},
    )
    (arguments.output / "paired_difference_draws.npz").chmod(0o440)
    np.savez_compressed(
        arguments.output / "model_metric_draws.npz",
        **{
            f"{model_id}--{metric_id}": draws[model_id][metric_id].astype(np.float64)
            for model_id in reference.MODEL_IDS
            for metric_id in METRIC_IDS
        },
    )
    (arguments.output / "model_metric_draws.npz").chmod(0o440)

    receipt = {
        "schema_version": "masld-bench-gse267145-histology-paired-contrast-v1",
        "analysis_revision_id": ANALYSIS_REVISION_ID,
        PRESPEC_SHA256_FIELD: arguments.prespec_sha256,
        "sealed_prespec_artifacts_sha256": arguments.sealed_prespec_artifacts_sha256,
        "scores_artifacts_sha256": SCORES_ARTIFACTS_SHA256,
        "scorer_source_sha256": SCORER_SOURCE_SHA256,
        "reference_evaluator_sha256": REFERENCE_EVALUATOR_SHA256,
        "outcomes_artifacts_sha256": OUTCOMES_ARTIFACTS_SHA256,
        "folds_artifacts_sha256": FOLDS_ARTIFACTS_SHA256,
        "bootstrap_indices_sha256": BOOTSTRAP_INDICES_SHA256,
        "bootstrap_indices_regenerated_and_verified": True,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "original_scoring_resamples_reused": True,
        "biological_resampling_unit": "participant",
        "participants": EXPECTED_PARTICIPANTS,
        "models_loaded": len(reference.MODEL_IDS),
        "prediction_bundles_read": 55,
        "prediction_bundles_mutated": False,
        "models_refitted": 0,
        "frozen_point_estimates_reproduced": reproduced,
        "declared_contrasts": len(rows),
        "observed_contrasts": len(observed),
        "not_applicable_contrasts": len(not_applicable),
        "per_replicate_difference_vectors_deposited": len(vectors),
        "metric_engine": "scripts.score_gse267145_histology_production_v4._bootstrap_metrics",
        "paired_cluster_bootstrap_used": False,
        "p_values_calculated": False,
        "bh_adjustment_calculated": False,
        "model_ranking_performed": False,
        "champion_claim_allowed": False,
        "external_transfer_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "confirmatory_inference_allowed": False,
        "between_sex_claim_allowed": False,
        "claim_mode": "single_cohort_participant_held_development_only",
        "status": "passed_paired_development_contrasts_no_confirmatory_claim",
    }
    write_json_exclusive(arguments.output / "contrast_receipt.json", receipt, mode=0o440)
    freeze_tree(
        arguments.output,
        {
            "artifact_class": "gse267145_histology_paired_development_contrasts",
            "scores_artifacts_sha256": SCORES_ARTIFACTS_SHA256,
            "sealed_prespec_artifacts_sha256": arguments.sealed_prespec_artifacts_sha256,
            "bootstrap_indices_sha256": BOOTSTRAP_INDICES_SHA256,
            "declared_contrasts": len(rows),
            "observed_contrasts": len(observed),
            "not_applicable_contrasts": len(not_applicable),
            "frozen_point_estimates_reproduced": reproduced,
            "models_refitted": 0,
            "p_values_calculated": False,
            "bh_adjustment_calculated": False,
            "model_ranking_performed": False,
            "status": "passed_paired_development_contrasts",
        },
    )
    print(json.dumps(receipt, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
