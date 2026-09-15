#!/usr/bin/env python3
"""Apply exact Atlas study-held baseline fits to outcome-blind GSE296875 RNA."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

import h5py
import numpy as np
from scipy import sparse

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file
from scripts.fit_predict_cell_baselines_study_50000 import (
    _log_normalize,
    _softmax,
    donor_class_weights,
    weighted_knn_probabilities,
)


VIEW_ID = "gse296875_rna_cell_state_external_development_7500_v1"
SOURCE_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
MODEL_IDS = (
    "hvg_pca_nearest_centroid",
    "hvg_pca_knn",
    "hvg_pca_logistic",
    "hvg_pca_elastic_net",
    "hvg_pca_linear_svm",
)
SEEDS = (20260824, 20260825, 20260826)
FOLDS = tuple(range(5))
TARGET_SUM = 10_000.0
SOURCE_REPRODUCTION_TOLERANCE = 5.0e-6


class BaselineTransferError(ValueError):
    """Raised when an outcome-blind transfer requirement differs."""


def _decode(values: Any) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or ()), list(reader)


def _write_predictions(path: Path, row_ids: np.ndarray, values: np.ndarray) -> None:
    fields = ("row_id", "predicted_class", *(f"probability::{label}" for label in ROSTER))
    predicted = np.argmax(values, axis=1)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, row_id in enumerate(row_ids):
            writer.writerow(
                {
                    "row_id": row_id,
                    "predicted_class": ROSTER[int(predicted[index])],
                    **{
                        f"probability::{label}": format(float(values[index, class_index]), ".17g")
                        for class_index, label in enumerate(ROSTER)
                    },
                }
            )


def _probabilities(values: np.ndarray, rows: int) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if (
        result.shape != (rows, len(ROSTER))
        or not np.all(np.isfinite(result))
        or np.any(result < 0.0)
    ):
        raise BaselineTransferError("probability tensor differs")
    sums = result.sum(axis=1, keepdims=True)
    if np.any(sums <= 0.0):
        raise BaselineTransferError("probability sum differs")
    result /= sums
    if not np.allclose(result.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12):
        raise BaselineTransferError("normalized probability sum differs")
    return result


def _linear_probability(x: np.ndarray, coefficient: np.ndarray, intercept: np.ndarray) -> np.ndarray:
    return _softmax(np.asarray(x, dtype=np.float64) @ coefficient.T + intercept)


def _centroid_probability(x: np.ndarray, centroids: np.ndarray, temperature: float) -> np.ndarray:
    distances = np.sum(
        (x[:, np.newaxis, :] - centroids[np.newaxis, :, :]) ** 2,
        axis=2,
    )
    return _softmax(-distances / float(temperature))


def _load_query(features: Path) -> tuple[sparse.csr_matrix, np.ndarray, list[str]]:
    with h5py.File(features / "rna_query.h5", "r") as handle:
        if set(handle["obs"].keys()) != {
            "row_id", "rna_status", "rna_observed_mask", "atac_observed_mask"
        }:
            raise BaselineTransferError("query observation firewall differs")
        row_ids = np.asarray(_decode(handle["obs/row_id"][:]), dtype=str)
        gene_ids = _decode(handle["rna/ensembl_id"][:])
        group = handle["rna/counts_csr"]
        matrix = sparse.csr_matrix(
            (
                np.asarray(group["data"][:]),
                np.asarray(group["indices"][:], dtype=np.int64),
                np.asarray(group["indptr"][:], dtype=np.int64),
            ),
            shape=tuple(int(value) for value in group["shape"][:]),
        )
        rna_mask = np.asarray(handle["obs/rna_observed_mask"][:])
        atac_mask = np.asarray(handle["obs/atac_observed_mask"][:])
    if (
        matrix.shape != (7_500, 36_601)
        or len(set(row_ids)) != 7_500
        or len(set(gene_ids)) != 36_601
        or np.any(rna_mask != 1)
        or np.any(atac_mask != 0)
        or np.any(matrix.data < 0)
        or np.any(matrix.data != np.floor(matrix.data))
    ):
        raise BaselineTransferError("query RNA matrix or missingness differs")
    return matrix, row_ids, gene_ids


def _load_source_predictions(
    baseline: Path,
    source_row_ids: np.ndarray,
) -> dict[tuple[int, str], np.ndarray]:
    expected_fields = [
        "row_id", "donor_id", "dataset", "outer_fold", "predicted_class",
        *(f"probability::{label}" for label in ROSTER),
    ]
    result: dict[tuple[int, str], np.ndarray] = {}
    for seed in SEEDS:
        for model_id in MODEL_IDS:
            fields, rows = _read_tsv(baseline / "predictions" / f"{model_id}--seed-{seed}.tsv")
            if fields != expected_fields or [row["row_id"] for row in rows] != source_row_ids.tolist():
                raise BaselineTransferError("source prediction row contract differs")
            result[(seed, model_id)] = _probabilities(
                np.asarray(
                    [[float(row[f"probability::{label}"]) for label in ROSTER] for row in rows]
                ),
                len(source_row_ids),
            )
    return result


def _selected_feature_axis(
    path: Path,
    source_gene_ids: list[str],
) -> tuple[list[str], np.ndarray]:
    fields, rows = _read_tsv(path)
    if fields != [
        "feature_index", "ensembl_id", "training_mean", "training_variance",
        "normalized_dispersion",
    ] or len(rows) != 2_000:
        raise BaselineTransferError("selected-feature schema differs")
    indices = np.asarray([int(row["feature_index"]) for row in rows], dtype=np.int64)
    ids = [row["ensembl_id"] for row in rows]
    if len(set(ids)) != 2_000 or ids != [source_gene_ids[index] for index in indices]:
        raise BaselineTransferError("selected-feature identity differs")
    return ids, indices


def run(args: argparse.Namespace) -> dict[str, Any]:
    import anndata

    if args.output.exists():
        raise BaselineTransferError("refusing to overwrite baseline transfer")
    features = args.features.resolve(strict=True)
    source = args.source.resolve(strict=True)
    split = args.split.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    for root, expected in (
        (features, args.expected_features_sha256),
        (source, args.expected_source_sha256),
        (split, args.expected_split_sha256),
        (baseline, args.expected_baseline_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise BaselineTransferError(f"input ARTIFACTS SHA-256 differs: {root}")
        verify_frozen_tree(root)
    feature_meta = verify_frozen_tree(features)["metadata"]
    source_meta = verify_frozen_tree(source)["metadata"]
    split_meta = verify_frozen_tree(split)["metadata"]
    baseline_meta = verify_frozen_tree(baseline)["metadata"]
    if (
        feature_meta.get("artifact_class")
        != "gse296875_rna_cell_state_external_development_features"
        or feature_meta.get("view_id") != VIEW_ID
        or feature_meta.get("labels_present") is not False
        or feature_meta.get("donor_ids_present") is not False
        or source_meta.get("subset_id") != SOURCE_VIEW_ID
        or source_meta.get("row_count") != 50_000
        or split_meta.get("split_id") != "resource_atlas_study_outer_5fold_v1"
        or split_meta.get("target_labels_used_for_assignment") is not False
        or baseline_meta.get("artifact_class") != "cell_state_task_native_predictions"
        or baseline_meta.get("dataset_view_id") != SOURCE_VIEW_ID
        or baseline_meta.get("metrics_calculated") is not False
        or baseline_meta.get("sealed_outcomes_read") is not False
    ):
        raise BaselineTransferError("input authority differs")

    query_raw, query_row_ids, query_gene_ids = _load_query(features)
    query_normalized = _log_normalize(query_raw, TARGET_SUM).tocsr()
    query_gene_index = {gene_id: index for index, gene_id in enumerate(query_gene_ids)}

    adata = anndata.read_h5ad(source / "resource_atlas_frozen_screen_50000.h5ad")
    source_row_ids = np.asarray(list(map(str, adata.obs["row_id"])), dtype=str)
    source_donors = np.asarray(list(map(str, adata.obs["donor_id"])), dtype=str)
    source_labels = np.asarray(list(map(str, adata.obs["broad_label"])), dtype=str)
    source_gene_ids = list(map(str, adata.var["ensembl_id"]))
    fields, split_rows = _read_tsv(split / "row_outer_folds.tsv")
    if (
        len(source_row_ids) != 50_000
        or [row["row_id"] for row in split_rows] != source_row_ids.tolist()
        or set(source_labels) != set(ROSTER)
    ):
        raise BaselineTransferError("source row contract differs")
    source_outer = np.asarray([int(row["outer_fold"]) for row in split_rows], dtype=np.int8)
    source_targets = np.asarray([ROSTER.index(label) for label in source_labels], dtype=np.int64)
    source_normalized = _log_normalize(adata.X, TARGET_SUM).tocsr()
    expected_source = _load_source_predictions(baseline, source_row_ids)

    query_by_model: dict[str, list[np.ndarray]] = {model_id: [] for model_id in MODEL_IDS}
    audit_rows: list[dict[str, Any]] = []
    maximum_error = 0.0
    missing_signatures: set[str] = set()
    for seed in SEEDS:
        for fold in FOLDS:
            fold_root = baseline / "folds" / f"seed{seed}" / f"fold{fold}"
            selected_ids, source_indices = _selected_feature_axis(
                fold_root / "selected_features.tsv", source_gene_ids
            )
            pca_mean = np.load(fold_root / "pca_mean.npy", allow_pickle=False)
            components = np.load(fold_root / "pca_components.npy", allow_pickle=False)
            if pca_mean.shape != (2_000,) or components.shape != (50, 2_000):
                raise BaselineTransferError("PCA state differs")

            source_dense = source_normalized[:, source_indices].toarray().astype(np.float32, copy=False)
            source_projection = (source_dense - pca_mean) @ components.T
            del source_dense
            query_dense = np.broadcast_to(pca_mean, (7_500, 2_000)).copy()
            present_positions: list[int] = []
            present_query_indices: list[int] = []
            missing_ids: list[str] = []
            for position, gene_id in enumerate(selected_ids):
                query_index = query_gene_index.get(gene_id)
                if query_index is None:
                    missing_ids.append(gene_id)
                else:
                    present_positions.append(position)
                    present_query_indices.append(query_index)
            query_dense[:, present_positions] = query_normalized[:, present_query_indices].toarray()
            query_projection = (query_dense - pca_mean) @ components.T
            del query_dense
            if not np.all(np.isfinite(source_projection)) or not np.all(np.isfinite(query_projection)):
                raise BaselineTransferError("PCA projection is nonfinite")

            train = source_outer != fold
            held = source_outer == fold
            weights = np.asarray(
                donor_class_weights(
                    source_donors[train].tolist(), source_labels[train].tolist()
                ),
                dtype=np.float64,
            )
            classifier_mean = np.load(
                fold_root / "classifier_projection_mean.npy", allow_pickle=False
            )
            classifier_scale = np.load(
                fold_root / "classifier_projection_scale.npy", allow_pickle=False
            )
            if classifier_mean.shape != (50,) or classifier_scale.shape != (50,) or np.any(classifier_scale <= 0):
                raise BaselineTransferError("classifier standardization state differs")
            source_standard = (source_projection - classifier_mean) / classifier_scale
            query_standard = (query_projection - classifier_mean) / classifier_scale

            centroid_state = json.loads((fold_root / "nearest_centroid.json").read_text())
            centroids = np.load(fold_root / "nearest_centroid_centroids.npy", allow_pickle=False)
            source_probabilities: dict[str, np.ndarray] = {
                "hvg_pca_nearest_centroid": _centroid_probability(
                    source_projection[held], centroids, float(centroid_state["temperature"])
                ),
                "hvg_pca_logistic": _linear_probability(
                    source_standard[held],
                    np.load(fold_root / "hvg_pca_logistic_coef.npy", allow_pickle=False),
                    np.load(fold_root / "hvg_pca_logistic_intercept.npy", allow_pickle=False),
                ),
                "hvg_pca_elastic_net": _linear_probability(
                    source_standard[held],
                    np.load(fold_root / "hvg_pca_elastic_net_coef.npy", allow_pickle=False),
                    np.load(fold_root / "hvg_pca_elastic_net_intercept.npy", allow_pickle=False),
                ),
            }
            query_probabilities: dict[str, np.ndarray] = {
                "hvg_pca_nearest_centroid": _centroid_probability(
                    query_projection, centroids, float(centroid_state["temperature"])
                ),
                "hvg_pca_logistic": _linear_probability(
                    query_standard,
                    np.load(fold_root / "hvg_pca_logistic_coef.npy", allow_pickle=False),
                    np.load(fold_root / "hvg_pca_logistic_intercept.npy", allow_pickle=False),
                ),
                "hvg_pca_elastic_net": _linear_probability(
                    query_standard,
                    np.load(fold_root / "hvg_pca_elastic_net_coef.npy", allow_pickle=False),
                    np.load(fold_root / "hvg_pca_elastic_net_intercept.npy", allow_pickle=False),
                ),
            }
            order = np.argsort(source_row_ids[train], kind="stable")
            source_probabilities["hvg_pca_knn"] = weighted_knn_probabilities(
                source_projection[train][order], source_targets[train][order], weights[order],
                source_projection[held], classes=len(ROSTER), neighbors=15,
                n_jobs=max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))),
            )
            query_probabilities["hvg_pca_knn"] = weighted_knn_probabilities(
                source_projection[train][order], source_targets[train][order], weights[order],
                query_projection, classes=len(ROSTER), neighbors=15,
                n_jobs=max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))),
            )
            svm_coef = np.load(fold_root / "svm_coef.npy", allow_pickle=False)
            svm_intercept = np.load(fold_root / "svm_intercept.npy", allow_pickle=False)
            calibration_coef = np.load(fold_root / "svm_calibration_coef.npy", allow_pickle=False)
            calibration_intercept = np.load(fold_root / "svm_calibration_intercept.npy", allow_pickle=False)
            source_scores = source_standard[held] @ svm_coef.T + svm_intercept
            query_scores = query_standard @ svm_coef.T + svm_intercept
            source_probabilities["hvg_pca_linear_svm"] = _linear_probability(
                source_scores, calibration_coef, calibration_intercept
            )
            query_probabilities["hvg_pca_linear_svm"] = _linear_probability(
                query_scores, calibration_coef, calibration_intercept
            )

            fold_errors: dict[str, float] = {}
            for model_id in MODEL_IDS:
                observed = expected_source[(seed, model_id)][held]
                reproduced = _probabilities(source_probabilities[model_id], int(held.sum()))
                error = float(np.max(np.abs(observed - reproduced)))
                if error > SOURCE_REPRODUCTION_TOLERANCE:
                    raise BaselineTransferError(
                        f"source prediction reproduction failed for {model_id}, seed {seed}, fold {fold}: {error}"
                    )
                fold_errors[model_id] = error
                maximum_error = max(maximum_error, error)
                query_by_model[model_id].append(
                    _probabilities(query_probabilities[model_id], 7_500)
                )
            missing_sha = canonical_sha256(missing_ids)
            missing_signatures.add(missing_sha)
            audit_rows.append(
                {
                    "seed": seed,
                    "outer_fold": fold,
                    "training_rows": int(train.sum()),
                    "training_donors": len(set(source_donors[train])),
                    "selected_features": 2_000,
                    "query_observed_selected_features": len(present_positions),
                    "query_structurally_missing_selected_features": len(missing_ids),
                    "query_structurally_missing_ids_sha256": missing_sha,
                    **{
                        f"source_reproduction_max_abs_error::{model_id}": format(fold_errors[model_id], ".17g")
                        for model_id in MODEL_IDS
                    },
                }
            )

    args.output.mkdir(mode=0o750)
    (args.output / "predictions").mkdir()
    for model_id, fits in query_by_model.items():
        if len(fits) != 15:
            raise BaselineTransferError("source-fit census differs")
        ensemble = _probabilities(np.mean(np.stack(fits), axis=0), 7_500)
        _write_predictions(args.output / "predictions" / f"{model_id}.tsv", query_row_ids, ensemble)
    with (args.output / "source_fit_audit.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(audit_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(audit_rows)
    receipt = {
        "schema_version": "masld-bench-gse296875-task-native-outer-ensemble-v1",
        "status": "pass_outcome_blind_predictions",
        "view_id": VIEW_ID,
        "models": list(MODEL_IDS),
        "rows": 7_500,
        "source_rows": 50_000,
        "source_fit_census_per_model": 15,
        "source_donors": 102,
        "source_studies": 7,
        "ensemble": "unweighted_probability_mean_over_three_seeds_by_five_study_outer_fits",
        "maximum_source_prediction_reproduction_error": maximum_error,
        "query_row_ids_sha256": canonical_sha256(query_row_ids.tolist()),
        "query_structurally_missing_feature_signatures": sorted(missing_signatures),
        "missing_feature_rule": "training_pca_mean_imputation_zero_contribution_after_centering_with_explicit_mask_audit",
        "observed_zero_rule": "retained_as_observed_zero_before_training_native_library_size_10000_log1p",
        "query_labels_read": False,
        "query_donor_ids_read": False,
        "query_atac_read": False,
        "atac_state": "structurally_missing",
        "target_adaptation_performed": False,
        "all_source_refit_performed": False,
        "metrics_calculated": False,
        "project_exposed_development_only": True,
        "external_or_champion_claim_allowed": False,
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(args.output / "prediction_receipt.json", receipt)
    freeze_tree(
        args.output,
        {
            "artifact_class": "gse296875_cell_state_outcome_blind_prediction_branch",
            "branch_id": "task_native_outer_ensemble",
            "view_id": VIEW_ID,
            "models": list(MODEL_IDS),
            "rows": 7_500,
            "source_fit_census_per_model": 15,
            "query_labels_read": False,
            "query_donor_ids_read": False,
            "query_atac_read": False,
            "target_adaptation_performed": False,
            "metrics_calculated": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--expected-features-sha256", required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--split", type=Path, required=True)
    parser.add_argument("--expected-split-sha256", required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--expected-baseline-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), sort_keys=True))
