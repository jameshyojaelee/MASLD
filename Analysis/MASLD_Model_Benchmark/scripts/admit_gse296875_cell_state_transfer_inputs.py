#!/usr/bin/env python3
"""Freeze outcome-blind authorities for GSE296875 cell-state transfer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


VIEW_ID = "gse296875_rna_cell_state_external_development_7500_v1"
SOURCE_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
BASELINE_MODELS = (
    "hvg_pca_nearest_centroid",
    "hvg_pca_knn",
    "hvg_pca_logistic",
    "hvg_pca_elastic_net",
    "hvg_pca_linear_svm",
)
SEEDS = (20260824, 20260825, 20260826)
FOLDS = tuple(range(5))
TF_HEADS = ("linear", "two_layer_mlp")


class TransferAdmissionError(ValueError):
    """Raised when an input authority is incomplete or outcome-exposed."""


def _metadata(root: Path, expected_sha256: str | None = None) -> dict[str, Any]:
    if expected_sha256 is not None and sha256_file(root / "ARTIFACTS.json") != expected_sha256:
        raise TransferAdmissionError(f"ARTIFACTS SHA-256 differs: {root}")
    return dict(verify_frozen_tree(root)["metadata"])


def admit(args: argparse.Namespace) -> dict[str, Any]:
    if args.output.exists():
        raise TransferAdmissionError("refusing to overwrite transfer admission")

    features = args.features.resolve(strict=True)
    verifier = args.verifier.resolve(strict=True)
    source = args.source.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    tf_embeddings = args.tf_embeddings.resolve(strict=True)
    tf_heads = args.tf_heads.resolve(strict=True)

    feature_meta = _metadata(features, args.expected_features_sha256)
    verifier_meta = _metadata(verifier, args.expected_verifier_sha256)
    source_meta = _metadata(source, args.expected_source_sha256)
    baseline_meta = _metadata(baseline)
    tf_embedding_meta = _metadata(tf_embeddings, args.expected_tf_embeddings_sha256)
    tf_head_meta = _metadata(tf_heads, args.expected_tf_heads_sha256)

    if (
        feature_meta.get("artifact_class")
        != "gse296875_rna_cell_state_external_development_features"
        or feature_meta.get("view_id") != VIEW_ID
        or feature_meta.get("rows") != 7_500
        or feature_meta.get("labels_present") is not False
        or feature_meta.get("donor_ids_present") is not False
        or feature_meta.get("sealed_outcomes_accessed") is not False
        or verifier_meta.get("artifact_class")
        != "gse296875_rna_cell_state_external_development_verification_campaign"
        or verifier_meta.get("sealed_outcomes_accessed") is not False
    ):
        raise TransferAdmissionError("GSE296875 feature or verifier firewall differs")
    if (
        source_meta.get("subset_id") != SOURCE_VIEW_ID
        or source_meta.get("row_count") != 50_000
        or source_meta.get("histology_read") is not False
        or source_meta.get("sealed_outcomes_read") is not False
        or baseline_meta.get("artifact_class") != "cell_state_task_native_predictions"
        or baseline_meta.get("dataset_view_id") != SOURCE_VIEW_ID
        or baseline_meta.get("rows") != 50_000
        or baseline_meta.get("donors") != 102
        or baseline_meta.get("studies") != 7
        or baseline_meta.get("metrics_calculated") is not False
        or baseline_meta.get("sealed_outcomes_read") is not False
    ):
        raise TransferAdmissionError("Atlas baseline authority differs")
    if (
        tf_embedding_meta.get("artifact_class")
        != "cell_foundation_frozen_screen_embeddings"
        or tf_embedding_meta.get("model_id") != "transcriptformer_tf_sapiens"
        or tf_embedding_meta.get("development_rows") != 50_000
        or tf_embedding_meta.get("sealed_outcomes_read") is not False
        or tf_head_meta.get("artifact_class")
        != "cell_foundation_common_head_shard_bundle"
        or tf_head_meta.get("logical_shards") != 30
        or tf_head_meta.get("head_ids") != list(TF_HEADS)
        or tf_head_meta.get("metrics_calculated") is not False
        or tf_head_meta.get("sealed_outcomes_read") is not False
    ):
        raise TransferAdmissionError("TF-Sapiens source-fit authority differs")

    baseline_files: list[dict[str, Any]] = []
    for seed in SEEDS:
        for fold in FOLDS:
            fold_root = baseline / "folds" / f"seed{seed}" / f"fold{fold}"
            required = [
                "fold_receipt.json",
                "selected_features.tsv",
                "pca_components.npy",
                "pca_mean.npy",
                "classifier_projection_mean.npy",
                "classifier_projection_scale.npy",
                "nearest_centroid_centroids.npy",
                "nearest_centroid.json",
                "hvg_pca_logistic_coef.npy",
                "hvg_pca_logistic_intercept.npy",
                "hvg_pca_elastic_net_coef.npy",
                "hvg_pca_elastic_net_intercept.npy",
                "svm_coef.npy",
                "svm_intercept.npy",
                "svm_calibration_coef.npy",
                "svm_calibration_intercept.npy",
            ]
            for relative in required:
                path = fold_root / relative
                if not path.is_file():
                    raise TransferAdmissionError(f"baseline source-fit file absent: {path}")
                baseline_files.append(
                    {"path": path.relative_to(baseline).as_posix(), "sha256": sha256_file(path)}
                )
    for model_id in BASELINE_MODELS:
        path = baseline / "predictions" / f"{model_id}.tsv"
        if not path.is_file():
            raise TransferAdmissionError(f"baseline prediction audit absent: {model_id}")

    tf_files: list[dict[str, Any]] = []
    for head_id in TF_HEADS:
        for seed in (1103, 1201, 1301):
            for fold in FOLDS:
                shard = tf_heads / "shards" / f"{head_id}__seed{seed}__fold{fold}"
                receipt = json.loads((shard / "prediction_receipt.json").read_text())
                if (
                    receipt.get("model_id") != "transcriptformer_tf_sapiens"
                    or receipt.get("head_id") != head_id
                    or receipt.get("screen_seed") != seed
                    or receipt.get("outer_fold") != fold
                    or receipt.get("outer_prediction_uses_final_refit") is not True
                    or receipt.get("input_artifacts_sha256", {}).get("split")
                    != "10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680"
                    or receipt.get("activation_fold_annotations_used_for_head_fitting") is not False
                    or receipt.get("sealed_outcomes_read") is not False
                ):
                    raise TransferAdmissionError("TF-Sapiens source-fit receipt differs")
                for relative in ("prediction_receipt.json", "predictions.tsv", "states/final_outer_training_refit.npz"):
                    path = shard / relative
                    tf_files.append(
                        {"path": path.relative_to(tf_heads).as_posix(), "sha256": sha256_file(path)}
                    )

    args.output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-gse296875-cell-state-transfer-input-authority-v1",
        "status": "pass_outcome_blind_transfer_inputs",
        "view_id": VIEW_ID,
        "source_view_id": SOURCE_VIEW_ID,
        "query_rows": 7_500,
        "source_rows": 50_000,
        "source_donors": 102,
        "source_studies": 7,
        "baseline_models": list(BASELINE_MODELS),
        "baseline_source_fit_census": len(SEEDS) * len(FOLDS),
        "tf_sapiens_heads": list(TF_HEADS),
        "tf_sapiens_source_fit_census_per_head": 15,
        "tf_sapiens_source_fit_census_total": 30,
        "common_transfer_rule": "unweighted_probability_mean_over_same_five_outer_folds_and_three_seeds",
        "all_source_refit_performed": False,
        "target_adaptation_performed": False,
        "query_evaluator_path_resolved": False,
        "query_labels_read": False,
        "query_donor_ids_read": False,
        "query_atac_read": False,
        "atac_missingness": "structurally_missing",
        "donor_fingerprint_overlap_status": "unresolved",
        "project_exposed_development_only": True,
        "external_or_champion_claim_allowed": False,
        "sealed_outcomes_read": False,
        "inputs": {
            "features": {"path": features.as_posix(), "artifacts_sha256": args.expected_features_sha256},
            "verifier": {"path": verifier.as_posix(), "artifacts_sha256": args.expected_verifier_sha256},
            "source": {"path": source.as_posix(), "artifacts_sha256": args.expected_source_sha256},
            "baseline": {"path": baseline.as_posix(), "artifacts_sha256": sha256_file(baseline / "ARTIFACTS.json")},
            "tf_embeddings": {"path": tf_embeddings.as_posix(), "artifacts_sha256": args.expected_tf_embeddings_sha256},
            "tf_heads": {"path": tf_heads.as_posix(), "artifacts_sha256": args.expected_tf_heads_sha256},
        },
        "baseline_source_fit_files": baseline_files,
        "tf_sapiens_source_fit_files": tf_files,
    }
    write_json_exclusive(args.output / "input_authority.json", receipt)
    freeze_tree(
        args.output,
        {
            "artifact_class": "gse296875_cell_state_transfer_input_authority",
            "view_id": VIEW_ID,
            "source_view_id": SOURCE_VIEW_ID,
            "query_rows": 7_500,
            "source_donors": 102,
            "baseline_source_fit_census": 15,
            "tf_sapiens_heads": list(TF_HEADS),
            "tf_sapiens_source_fit_census_per_head": 15,
            "tf_sapiens_source_fit_census_total": 30,
            "query_labels_read": False,
            "query_donor_ids_read": False,
            "query_atac_read": False,
            "external_or_champion_claim_allowed": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--expected-features-sha256", required=True)
    parser.add_argument("--verifier", type=Path, required=True)
    parser.add_argument("--expected-verifier-sha256", required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--tf-embeddings", type=Path, required=True)
    parser.add_argument("--expected-tf-embeddings-sha256", required=True)
    parser.add_argument("--tf-heads", type=Path, required=True)
    parser.add_argument("--expected-tf-heads-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(admit(parse_args()), sort_keys=True))
