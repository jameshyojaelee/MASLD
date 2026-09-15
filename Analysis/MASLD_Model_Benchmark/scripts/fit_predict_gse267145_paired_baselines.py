#!/usr/bin/env python3
"""Fit prediction-only GSE267145 RNA-to-H3 baselines for one outer fold."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.utils.extmath import randomized_svd

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
EXPECTED_RNA_FEATURES = 42_163
EXPECTED_H3_FEATURES = 96_460
EXPECTED_FOLD_COUNTS = {0: 21, 1: 21, 2: 21, 3: 19, 4: 17}
RNA_VARIANCE_FEATURES = 4_096
RNA_PCS = 20
H3_PCS = 30
RIDGE_ALPHA = 10.0
RRR_ALPHA = 1.0
RRR_RANK = 10
PLS_COMPONENTS = 10
SCCA_RNA_FEATURES = 512
SCCA_H3_FEATURES = 2_048
SCCA_COMPONENTS = 5
SCCA_RNA_NONZERO = 64
SCCA_H3_NONZERO = 128
PROFILE_FLOOR = 1.0e-12
LIBRARY_SCALE = 1.0e6


class PairedBaselineFitError(ValueError):
    """Raised when a fit would not meet the outer-fold or numeric requirements."""


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PairedBaselineFitError(f"TSV lacks a header: {path}")
        return list(reader.fieldnames), list(reader)


def library_log1p(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or not np.all(np.isfinite(array)) or np.any(array < 0.0):
        raise PairedBaselineFitError("RNA values must be finite nonnegative matrices")
    totals = array.sum(axis=1)
    if np.any(totals <= 0.0):
        raise PairedBaselineFitError("RNA participant library total is nonpositive")
    return np.log1p(array / totals[:, np.newaxis] * LIBRARY_SCALE)


def h3_hellinger(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    array = np.asarray(counts)
    if array.ndim != 2 or not np.issubdtype(array.dtype, np.integer) or np.any(array < 0):
        raise PairedBaselineFitError("H3 values must be nonnegative integer matrices")
    totals = array.sum(axis=1, dtype=np.float64)
    if np.any(totals <= 0.0):
        raise PairedBaselineFitError("H3 participant library total is nonpositive")
    profiles = array.astype(np.float64) / totals[:, np.newaxis]
    return np.sqrt(profiles), totals


def normalize_profile(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or not np.all(np.isfinite(array)):
        raise PairedBaselineFitError("predicted H3 profile matrix is invalid")
    array = np.maximum(array, 0.0) ** 2
    array += PROFILE_FLOOR
    totals = array.sum(axis=1)
    if np.any(totals <= 0.0):
        raise PairedBaselineFitError("predicted H3 profile total is nonpositive")
    result = array / totals[:, np.newaxis]
    if not np.allclose(result.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12):
        raise PairedBaselineFitError("predicted H3 profiles do not sum to one")
    return result


def top_variance_indices(values: np.ndarray, maximum: int) -> np.ndarray:
    variance = np.var(values, axis=0, ddof=1)
    eligible = np.flatnonzero(np.isfinite(variance) & (variance > 0.0))
    if not len(eligible):
        raise PairedBaselineFitError("no positive-variance features remain")
    take = min(maximum, len(eligible))
    order = np.lexsort((eligible, -variance[eligible]))[:take]
    return eligible[order].astype(np.int32)


def standardized_columns(
    training: np.ndarray, query: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mean = np.mean(training, axis=0)
    scale = np.std(training, axis=0, ddof=1)
    scale = np.where(scale > 0.0, scale, 1.0)
    return (
        (training - mean) / scale,
        (query - mean) / scale,
        mean,
        scale,
    )


def rna_pca_scores(
    training: np.ndarray, query: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = top_variance_indices(training, RNA_VARIANCE_FEATURES)
    train_scaled, query_scaled, _, _ = standardized_columns(
        training[:, selected], query[:, selected]
    )
    components = min(RNA_PCS, len(training) - 1, train_scaled.shape[1])
    if components < 2:
        raise PairedBaselineFitError("RNA PCA rank is insufficient")
    _, _, right = np.linalg.svd(train_scaled, full_matrices=False)
    loadings = right[:components]
    return train_scaled @ loadings.T, query_scaled @ loadings.T, selected


def target_pca(
    training_hellinger: np.ndarray, *, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.mean(training_hellinger, axis=0)
    centered = training_hellinger - mean
    components = min(H3_PCS, len(training_hellinger) - 1, centered.shape[1])
    left, singular, right = randomized_svd(
        centered,
        n_components=components,
        n_iter=7,
        random_state=seed,
        flip_sign=True,
    )
    return left * singular, right, mean


def ridge_predict(
    train_x: np.ndarray,
    query_x: np.ndarray,
    train_y: np.ndarray,
    *,
    alpha: float,
) -> np.ndarray:
    x_mean = np.mean(train_x, axis=0)
    y_mean = np.mean(train_y, axis=0)
    centered_x = train_x - x_mean
    coefficient = np.linalg.solve(
        centered_x.T @ centered_x + alpha * np.eye(centered_x.shape[1]),
        centered_x.T @ (train_y - y_mean),
    )
    return (query_x - x_mean) @ coefficient + y_mean


def reduced_rank_predict(
    train_x: np.ndarray,
    query_x: np.ndarray,
    train_y: np.ndarray,
) -> np.ndarray:
    x_mean = np.mean(train_x, axis=0)
    y_mean = np.mean(train_y, axis=0)
    centered_x = train_x - x_mean
    centered_y = train_y - y_mean
    coefficient = np.linalg.solve(
        centered_x.T @ centered_x + RRR_ALPHA * np.eye(centered_x.shape[1]),
        centered_x.T @ centered_y,
    )
    fitted = centered_x @ coefficient
    rank = min(RRR_RANK, len(train_x) - 1, train_x.shape[1])
    _, _, right = randomized_svd(
        fitted,
        n_components=rank,
        n_iter=7,
        random_state=267145,
        flip_sign=True,
    )
    return ((query_x - x_mean) @ coefficient @ right.T) @ right + y_mean


def _soft_sparse(vector: np.ndarray, nonzero: int) -> np.ndarray:
    values = np.asarray(vector, dtype=np.float64)
    if nonzero >= len(values):
        result = values.copy()
    else:
        threshold = np.partition(np.abs(values), -nonzero - 1)[-nonzero - 1]
        result = np.sign(values) * np.maximum(np.abs(values) - threshold, 0.0)
    norm = np.linalg.norm(result)
    if not np.isfinite(norm) or norm <= 0.0:
        raise PairedBaselineFitError("sparse canonical loading collapsed")
    return result / norm


def sparse_cca_loadings(
    cross_covariance: np.ndarray,
    *,
    components: int = SCCA_COMPONENTS,
    x_nonzero: int = SCCA_RNA_NONZERO,
    y_nonzero: int = SCCA_H3_NONZERO,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    residual = np.asarray(cross_covariance, dtype=np.float64).copy()
    x_loadings: list[np.ndarray] = []
    y_loadings: list[np.ndarray] = []
    iterations: list[int] = []
    for _ in range(components):
        left, _, right = randomized_svd(
            residual,
            n_components=1,
            n_iter=7,
            random_state=267145 + len(x_loadings),
            flip_sign=True,
        )
        u = _soft_sparse(left[:, 0], min(x_nonzero, residual.shape[0]))
        v = _soft_sparse(right[0], min(y_nonzero, residual.shape[1]))
        for iteration in range(1, 101):
            next_u = _soft_sparse(
                residual @ v, min(x_nonzero, residual.shape[0])
            )
            next_v = _soft_sparse(
                residual.T @ next_u, min(y_nonzero, residual.shape[1])
            )
            delta = max(np.linalg.norm(next_u - u), np.linalg.norm(next_v - v))
            u, v = next_u, next_v
            if delta <= 1.0e-8:
                break
        covariance = float(u @ residual @ v)
        if not np.isfinite(covariance) or covariance <= 0.0:
            raise PairedBaselineFitError("sparse canonical covariance is nonpositive")
        x_loadings.append(u)
        y_loadings.append(v)
        iterations.append(iteration)
        residual -= covariance * np.outer(u, v)
    return np.column_stack(x_loadings), np.column_stack(y_loadings), iterations


def fit_fold(*, view: Path, output: Path, expected_view_sha256: str) -> dict[str, Any]:
    if output.exists():
        raise PairedBaselineFitError("refusing to overwrite fold predictions")
    if sha256_file(view / "ARTIFACTS.json") != expected_view_sha256:
        raise PairedBaselineFitError("model view ARTIFACTS SHA-256 differs")
    metadata = verify_frozen_tree(view)["metadata"]
    receipt = json.loads((view / "receipt.json").read_text(encoding="utf-8"))
    outer_fold = int(receipt.get("outer_fold", -1))
    if (
        metadata.get("artifact_class") != "gse267145_paired_model_fold_view"
        or metadata.get("outer_fold") != outer_fold
        or outer_fold not in range(5)
        or receipt.get("query_h3k27ac_included") is not False
        or receipt.get("held_h3_library_totals_included") is not False
        or receipt.get("outcome_or_participant_covariates_included") is not False
        or receipt.get("masked_retired_rna_genes") != 1_122
        or receipt.get("h3_feature_identity") != "opaque_source_feature_key"
        or receipt.get("coordinate_semantics_inferred") is not False
    ):
        raise PairedBaselineFitError("model view firewall metadata differs")
    if any("query_h3" in path.name.lower() for path in view.iterdir()):
        raise PairedBaselineFitError("model view contains held H3 material")
    train_fields, training_rows = read_tsv(view / "training_participants.tsv")
    query_fields, query_rows = read_tsv(view / "query_participants.tsv")
    expected_fields = ["participant_index", "participant_id", "outer_fold"]
    if train_fields != expected_fields or query_fields != expected_fields:
        raise PairedBaselineFitError("participant table schema differs")
    if (
        len(training_rows) + len(query_rows) != 99
        or len(query_rows) != EXPECTED_FOLD_COUNTS[outer_fold]
        or {row["participant_id"] for row in training_rows}
        & {row["participant_id"] for row in query_rows}
        or any(int(row["outer_fold"]) == outer_fold for row in training_rows)
        or any(int(row["outer_fold"]) != outer_fold for row in query_rows)
    ):
        raise PairedBaselineFitError("outer participant partition differs")

    train_rna = np.load(view / "training_rna_values.npy", allow_pickle=False)
    query_rna = np.load(view / "query_rna_values.npy", allow_pickle=False)
    train_h3 = np.load(view / "training_h3k27ac_counts.npy", allow_pickle=False)
    if (
        train_rna.shape != (len(training_rows), EXPECTED_RNA_FEATURES)
        or query_rna.shape != (len(query_rows), EXPECTED_RNA_FEATURES)
        or train_h3.shape != (len(training_rows), EXPECTED_H3_FEATURES)
        or train_h3.dtype != np.uint32
    ):
        raise PairedBaselineFitError("model view array shape differs")
    transformed_train_rna = library_log1p(train_rna)
    transformed_query_rna = library_log1p(query_rna)
    train_h3_sqrt, _ = h3_hellinger(train_h3)
    train_h3_profiles = train_h3_sqrt**2
    x_train, x_query, rna_pca_indices = rna_pca_scores(
        transformed_train_rna, transformed_query_rna
    )
    target_scores, target_loadings, target_mean = target_pca(
        train_h3_sqrt, seed=267145 + outer_fold
    )

    predictions: dict[str, np.ndarray] = {}
    mean_profile = np.mean(train_h3_profiles, axis=0)
    mean_profile /= mean_profile.sum()
    predictions["training_mean_h3_profile"] = np.repeat(
        normalize_profile(np.sqrt(mean_profile)[np.newaxis, :]),
        len(query_rows),
        axis=0,
    )
    pca_scores = ridge_predict(
        x_train, x_query, target_scores, alpha=RIDGE_ALPHA
    )
    predictions["pca_ridge"] = normalize_profile(
        target_mean + pca_scores @ target_loadings
    )
    predictions["reduced_rank_regression"] = normalize_profile(
        reduced_rank_predict(x_train, x_query, train_h3_sqrt)
    )
    pls_components = min(
        PLS_COMPONENTS,
        x_train.shape[1],
        len(x_train) - 1,
    )
    pls = PLSRegression(
        n_components=pls_components,
        scale=False,
        max_iter=500,
        tol=1.0e-7,
        copy=True,
    )
    pls.fit(x_train, train_h3_sqrt)
    if any(iterations >= 500 for iterations in pls.n_iter_):
        raise PairedBaselineFitError("PLS2 reached the fixed iteration limit")
    predictions["pls2"] = normalize_profile(pls.predict(x_query))

    scca_rna_indices = top_variance_indices(
        transformed_train_rna, SCCA_RNA_FEATURES
    )
    scca_h3_indices = top_variance_indices(train_h3_sqrt, SCCA_H3_FEATURES)
    scca_x_train, scca_x_query, _, _ = standardized_columns(
        transformed_train_rna[:, scca_rna_indices],
        transformed_query_rna[:, scca_rna_indices],
    )
    scca_y_train, _, scca_h3_mean, scca_h3_scale = standardized_columns(
        train_h3_sqrt[:, scca_h3_indices],
        train_h3_sqrt[:, scca_h3_indices],
    )
    x_weights, h3_weights, scca_iterations = sparse_cca_loadings(
        scca_x_train.T @ scca_y_train / (len(scca_x_train) - 1)
    )
    query_embeddings = scca_x_query @ x_weights
    if not np.all(np.isfinite(query_embeddings)):
        raise PairedBaselineFitError("sparse CCA query embeddings are invalid")

    output.mkdir(mode=0o750)
    prediction_root = output / "predictions"
    state_root = output / "sparse_cca"
    prediction_root.mkdir()
    state_root.mkdir()
    for model_id, values in predictions.items():
        if (
            values.shape != (len(query_rows), EXPECTED_H3_FEATURES)
            or not np.all(np.isfinite(values))
            or np.any(values <= 0.0)
            or not np.allclose(values.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-7)
        ):
            raise PairedBaselineFitError(f"{model_id} profile predictions differ")
        np.save(
            prediction_root / f"{model_id}.npy",
            values.astype(np.float32),
            allow_pickle=False,
        )
    np.save(
        state_root / "query_embeddings.npy",
        query_embeddings.astype(np.float32),
        allow_pickle=False,
    )
    with (state_root / "candidate_transform.npz").open("xb") as handle:
        np.savez_compressed(
            handle,
            h3_feature_indices=scca_h3_indices.astype(np.int32),
            h3_mean=scca_h3_mean.astype(np.float32),
            h3_scale=scca_h3_scale.astype(np.float32),
            h3_weights=h3_weights.astype(np.float32),
        )
    write_tsv_fields = ["participant_index", "participant_id", "outer_fold"]
    with (output / "query_participants.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=write_tsv_fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(query_rows)
    result = {
        "schema_version": "masld-bench-gse267145-paired-baseline-fold-v1",
        "status": "pass_development_predictions",
        "campaign_stage": "task_native_baseline_screen",
        "outer_fold": outer_fold,
        "training_participants": len(training_rows),
        "query_participants": len(query_rows),
        "model_ids": list(MODEL_IDS),
        "profile_prediction_models": list(PROFILE_MODEL_IDS),
        "retrieval_models": list(MODEL_IDS),
        "sparse_cca_profile_prediction_applicable": False,
        "sparse_cca_profile_disposition": "not_applicable_symmetric_association_method_not_forced_into_directional_profile_regression",
        "sparse_cca_retrieval_applicable": True,
        "preprocessing": {
            "rna": "participant_local_library_scale_1e6_then_log1p",
            "rna_variance_features": RNA_VARIANCE_FEATURES,
            "rna_pcs": x_train.shape[1],
            "h3_target": "participant_count_composition_hellinger",
            "h3_pcs": target_scores.shape[1],
            "pls2_response_space": "full_h3_hellinger_profile",
            "profile_floor": PROFILE_FLOOR,
            "all_learned_parameters_fit_on_outer_training_only": True,
        },
        "fixed_hyperparameters": {
            "pca_ridge_alpha": RIDGE_ALPHA,
            "reduced_rank_alpha": RRR_ALPHA,
            "reduced_rank_rank": RRR_RANK,
            "pls2_components": pls_components,
            "pls2_iterations": [int(value) for value in pls.n_iter_],
            "sparse_cca_components": SCCA_COMPONENTS,
            "sparse_cca_rna_screen": len(scca_rna_indices),
            "sparse_cca_h3_screen": len(scca_h3_indices),
            "sparse_cca_rna_nonzero_per_component_max": SCCA_RNA_NONZERO,
            "sparse_cca_h3_nonzero_per_component_max": SCCA_H3_NONZERO,
            "sparse_cca_iterations": scca_iterations,
        },
        "fixed_hyperparameters_selected_without_source_outcomes": True,
        "random_seeds": {
            "target_pca": 267145 + outer_fold,
            "reduced_rank_svd": 267145,
            "sparse_cca_component_base": 267145,
        },
        "rna_pca_feature_indices_sha256": canonical_sha256(rna_pca_indices.tolist()),
        "sparse_cca_rna_feature_indices_sha256": canonical_sha256(
            scca_rna_indices.tolist()
        ),
        "sparse_cca_h3_opaque_feature_indices_sha256": canonical_sha256(
            scca_h3_indices.tolist()
        ),
        "h3_feature_identity": "opaque_source_feature_index",
        "coordinate_interpretation_used": False,
        "held_participant_h3_read": False,
        "held_participant_h3_library_totals_read": False,
        "histology_sex_age_nas_fibrosis_read": False,
        "models_fit": True,
        "metrics_calculated": False,
        "prediction_tables_contain_observed_h3": False,
        "artifact_release_state": "internal_only_pending_model_specific_rights_review",
        "derivative_weight_release_requires_review": True,
        "champion_claim_allowed": False,
        "input_view_artifacts_sha256": expected_view_sha256,
    }
    write_json_exclusive(output / "prediction_receipt.json", result)
    freeze_tree(
        output,
        {
            "artifact_class": "gse267145_paired_baseline_fold_predictions",
            "outer_fold": outer_fold,
            "models": len(MODEL_IDS),
            "query_participants": len(query_rows),
            "held_participant_h3_read": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    return result


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--view", required=True, type=Path)
    value.add_argument("--view-artifacts-sha256", required=True)
    value.add_argument("--output", required=True, type=Path)
    return value


def main() -> int:
    arguments = parser().parse_args()
    result = fit_fold(
        view=arguments.view,
        output=arguments.output,
        expected_view_sha256=arguments.view_artifacts_sha256,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
