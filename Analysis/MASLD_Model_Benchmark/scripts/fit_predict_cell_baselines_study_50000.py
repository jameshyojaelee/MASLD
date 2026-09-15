#!/usr/bin/env python3
"""Fit study-held-out classical cell-state baselines and freeze predictions.

This scientific adapter reads development labels to fit models but never scores
them. Metrics remain the responsibility of an independent evaluator.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence
import warnings

import numpy as np

ADAPTER_ROOT = Path(__file__).parents[1] / "src/masld_bench/adapters"
sys.path.insert(0, ADAPTER_ROOT.as_posix())
from hvg_pca_logistic import (  # type: ignore[import-not-found]  # noqa: E402
    _log_normalize,
    _select_hvgs,
    _softmax,
    donor_class_weights,
)


DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
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
N_TOP_HVG = 2_000
N_PCA = 50
TARGET_SUM = 10_000.0
KNN_NEIGHBORS = 15
SEEDS = (20260824, 20260825, 20260826)
LOGISTIC_MAX_ITER = 2_000
ELASTIC_NET_MAX_ITER = 5_000
ELASTIC_NET_TOL = 1.0e-4
ELASTIC_NET_N_ITER_NO_CHANGE = 20
ELASTIC_NET_ALPHAS = (1.0e-5, 1.0e-4, 1.0e-3)
ELASTIC_NET_L1_RATIOS = (0.15, 0.5, 0.85)
ELASTIC_NET_MIN_VALID_CANDIDATES = 2


class BaselineScreenError(ValueError):
    """Raised when a baseline input, fit, or prediction violates requirement."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def write_json_exclusive(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(
            value,
            handle,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        handle.write("\n")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise BaselineScreenError(f"TSV lacks a header: {path}")
        return list(reader)


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
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


def nearest_centroid_probabilities(
    train_x: np.ndarray,
    train_y: np.ndarray,
    train_weights: np.ndarray,
    test_x: np.ndarray,
    *,
    classes: int,
) -> tuple[np.ndarray, dict[str, np.ndarray | float]]:
    centroids = np.vstack(
        [
            np.average(train_x[train_y == class_id], axis=0, weights=train_weights[train_y == class_id])
            for class_id in range(classes)
        ]
    )
    within = np.sum((train_x - centroids[train_y]) ** 2, axis=1)
    temperature = max(float(np.average(within, weights=train_weights)), 1.0e-12)
    distances = np.sum(
        (test_x[:, np.newaxis, :] - centroids[np.newaxis, :, :]) ** 2,
        axis=2,
    )
    return _softmax(-distances / temperature), {
        "centroids": centroids,
        "temperature": temperature,
    }


def weighted_knn_probabilities(
    train_x: np.ndarray,
    train_y: np.ndarray,
    train_weights: np.ndarray,
    test_x: np.ndarray,
    *,
    classes: int,
    neighbors: int,
    n_jobs: int,
) -> np.ndarray:
    from sklearn.neighbors import NearestNeighbors

    if not 0 < neighbors <= len(train_x):
        raise BaselineScreenError("kNN neighbor count differs")
    search = NearestNeighbors(n_neighbors=neighbors, algorithm="auto", n_jobs=n_jobs)
    search.fit(train_x)
    _, indices = search.kneighbors(test_x, return_distance=True)
    result = np.zeros((len(test_x), classes), dtype=np.float64)
    for row_index, nearest in enumerate(indices):
        votes = np.bincount(
            train_y[nearest], weights=train_weights[nearest], minlength=classes
        )
        result[row_index] = votes / votes.sum()
    return result


def _validate_probabilities(values: np.ndarray, rows: int) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.shape != (rows, len(ROSTER)) or not np.all(np.isfinite(result)) or np.any(result < 0.0):
        raise BaselineScreenError("classifier probabilities differ")
    row_sums = result.sum(axis=1)
    if np.any(row_sums <= 0.0) or not np.allclose(
        row_sums, 1.0, rtol=1.0e-5, atol=1.0e-7
    ):
        raise BaselineScreenError("classifier probability sums differ")
    result = result / row_sums[:, np.newaxis]
    if not np.allclose(result.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12):
        raise BaselineScreenError("normalized classifier probability sums differ")
    return result


def _save_array(path: Path, value: Any) -> None:
    with path.open("xb") as handle:
        np.save(handle, np.asarray(value), allow_pickle=False)


def standardize_projection(
    train_x: np.ndarray,
    test_x: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit an unweighted z-score transform on training projections only."""
    train = np.asarray(train_x, dtype=np.float64)
    test = np.asarray(test_x, dtype=np.float64)
    if (
        train.ndim != 2
        or test.ndim != 2
        or train.shape[1] != test.shape[1]
        or len(train) == 0
        or not np.all(np.isfinite(train))
        or not np.all(np.isfinite(test))
    ):
        raise BaselineScreenError("classifier projection standardization input differs")
    mean = train.mean(axis=0)
    scale = train.std(axis=0)
    scale[scale <= np.finfo(np.float64).eps] = 1.0
    standardized_train = (train - mean) / scale
    standardized_test = (test - mean) / scale
    if not np.all(np.isfinite(standardized_train)) or not np.all(
        np.isfinite(standardized_test)
    ):
        raise BaselineScreenError("classifier projection standardization is nonfinite")
    return standardized_train, standardized_test, mean, scale


def validate_integer_like_counts(matrix: Any, *, chunk_size: int = 1_000_000) -> None:
    values = matrix.data if hasattr(matrix, "data") else np.asarray(matrix).reshape(-1)
    if chunk_size < 1:
        raise BaselineScreenError("raw-count validation chunk size differs")
    for start in range(0, len(values), chunk_size):
        block = np.asarray(values[start : start + chunk_size])
        if (
            not np.all(np.isfinite(block))
            or np.any(block < 0)
            or np.any(block != np.floor(block))
        ):
            raise BaselineScreenError(
                "Atlas matrix is not nonnegative integer-valued raw UMI counts"
            )


def fit_sgd_elastic_net(
    train_x: np.ndarray,
    train_y: np.ndarray,
    sample_weight: np.ndarray,
    *,
    alpha: float,
    l1_ratio: float,
    seed: int,
) -> tuple[Any, dict[str, Any]]:
    """Fit one deterministic, weighted elastic-net log-loss candidate."""
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import SGDClassifier

    if (
        train_x.ndim != 2
        or len(train_x) != len(train_y)
        or len(train_y) != len(sample_weight)
        or set(train_y) != set(range(len(ROSTER)))
        or not np.all(np.isfinite(train_x))
        or not np.all(np.isfinite(sample_weight))
        or np.any(sample_weight <= 0.0)
        or alpha not in ELASTIC_NET_ALPHAS
        or l1_ratio not in ELASTIC_NET_L1_RATIOS
    ):
        raise BaselineScreenError("SGD elastic-net fitting contract differs")
    classifier = SGDClassifier(
        loss="log_loss",
        penalty="elasticnet",
        alpha=alpha,
        l1_ratio=l1_ratio,
        fit_intercept=True,
        max_iter=ELASTIC_NET_MAX_ITER,
        tol=ELASTIC_NET_TOL,
        shuffle=True,
        random_state=seed,
        learning_rate="optimal",
        early_stopping=False,
        n_iter_no_change=ELASTIC_NET_N_ITER_NO_CHANGE,
        average=True,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        classifier.fit(train_x, train_y, sample_weight=sample_weight)
    convergence_warnings = [
        str(item.message) for item in caught if issubclass(item.category, ConvergenceWarning)
    ]
    converged = (
        tuple(classifier.classes_) == tuple(range(len(ROSTER)))
        and int(classifier.n_iter_) < ELASTIC_NET_MAX_ITER
        and not convergence_warnings
        and np.all(np.isfinite(classifier.coef_))
        and np.all(np.isfinite(classifier.intercept_))
    )
    state = {
        "solver": "SGDClassifier",
        "loss": "log_loss",
        "penalty": "elasticnet",
        "alpha": alpha,
        "l1_ratio": l1_ratio,
        "iterations": int(classifier.n_iter_),
        "max_iter": ELASTIC_NET_MAX_ITER,
        "tolerance": ELASTIC_NET_TOL,
        "n_iter_no_change": ELASTIC_NET_N_ITER_NO_CHANGE,
        "learning_rate": "optimal",
        "average": True,
        "convergence_warnings": convergence_warnings,
        "converged_before_ceiling": bool(converged),
        "coefficient_zero_fraction": float(np.mean(classifier.coef_ == 0.0)),
    }
    return classifier, state


def prepare_inner_projections(
    *,
    raw_train: Any,
    gene_ids: Sequence[str],
    train_y: np.ndarray,
    train_donors: np.ndarray,
    inner_assignment: Mapping[str, int],
    fold: int,
    seed: int,
) -> list[dict[str, Any]]:
    """Fit all feature/PCA transforms inside each inner training partition."""
    from sklearn.decomposition import PCA

    prepared: list[dict[str, Any]] = []
    for inner_fold in range(5):
        validation = np.asarray(
            [inner_assignment[str(donor)] == inner_fold for donor in train_donors]
        )
        fitting = ~validation
        if (
            not validation.any()
            or set(train_y[fitting]) != set(range(len(ROSTER)))
            or set(train_y[validation]) != set(range(len(ROSTER)))
        ):
            raise BaselineScreenError("inner fold lacks rows or classes")
        normalized_fitting = _log_normalize(raw_train[fitting], TARGET_SUM)
        inner_selected, _ = _select_hvgs(
            normalized_fitting, gene_ids, n_top=N_TOP_HVG, n_bins=20
        )
        fitting_dense = normalized_fitting[:, inner_selected].toarray().astype(
            np.float32, copy=False
        )
        normalized_validation = _log_normalize(raw_train[validation], TARGET_SUM)
        validation_dense = normalized_validation[:, inner_selected].toarray().astype(
            np.float32, copy=False
        )
        inner_pca = PCA(
            n_components=N_PCA,
            svd_solver="randomized",
            whiten=False,
            random_state=seed + fold * 101 + inner_fold,
        )
        fitting_x = inner_pca.fit_transform(fitting_dense).astype(np.float32, copy=False)
        validation_x = inner_pca.transform(validation_dense).astype(np.float32, copy=False)
        fitting_x, validation_x, _, _ = standardize_projection(fitting_x, validation_x)
        fitting_weights = np.asarray(
            donor_class_weights(
                list(train_donors[fitting]),
                [ROSTER[index] for index in train_y[fitting]],
            ),
            dtype=np.float64,
        )
        validation_weights = np.asarray(
            donor_class_weights(
                list(train_donors[validation]),
                [ROSTER[index] for index in train_y[validation]],
            ),
            dtype=np.float64,
        )
        prepared.append(
            {
                "inner_fold": inner_fold,
                "fitting": fitting,
                "validation": validation,
                "fitting_x": fitting_x,
                "validation_x": validation_x,
                "fitting_weights": fitting_weights,
                "validation_weights": validation_weights,
                "selected_feature_ids_sha256": canonical_sha256(
                    [gene_ids[index] for index in inner_selected]
                ),
            }
        )
    return prepared


def select_sgd_elastic_net(
    prepared: Sequence[Mapping[str, Any]],
    train_y: np.ndarray,
    *,
    fold: int,
    seed: int,
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    """Select hyperparameters using only donor-grouped inner validation rows."""
    from sklearn.metrics import log_loss

    candidate_audit: list[dict[str, Any]] = []
    for candidate_index, (alpha, l1_ratio) in enumerate(
        (alpha, l1_ratio)
        for alpha in ELASTIC_NET_ALPHAS
        for l1_ratio in ELASTIC_NET_L1_RATIOS
    ):
        fold_losses: list[float] = []
        fit_states: list[dict[str, Any]] = []
        for item in prepared:
            fitting = np.asarray(item["fitting"], dtype=bool)
            validation = np.asarray(item["validation"], dtype=bool)
            classifier, state = fit_sgd_elastic_net(
                np.asarray(item["fitting_x"]),
                train_y[fitting],
                np.asarray(item["fitting_weights"]),
                alpha=alpha,
                l1_ratio=l1_ratio,
                seed=seed + fold * 10_000 + int(item["inner_fold"]) * 100 + candidate_index,
            )
            fit_states.append(state)
            if not state["converged_before_ceiling"]:
                continue
            probabilities = _validate_probabilities(
                classifier.predict_proba(np.asarray(item["validation_x"])),
                int(validation.sum()),
            )
            fold_losses.append(
                float(
                    log_loss(
                        train_y[validation],
                        probabilities,
                        labels=np.arange(len(ROSTER)),
                        sample_weight=np.asarray(item["validation_weights"]),
                    )
                )
            )
        valid = len(fold_losses) == len(prepared)
        candidate_audit.append(
            {
                "alpha": alpha,
                "l1_ratio": l1_ratio,
                "valid_all_inner_folds": valid,
                "inner_fold_log_losses": fold_losses,
                "mean_inner_log_loss": float(np.mean(fold_losses)) if valid else None,
                "fit_states": fit_states,
            }
        )
    valid_candidates = [item for item in candidate_audit if item["valid_all_inner_folds"]]
    if len(valid_candidates) < ELASTIC_NET_MIN_VALID_CANDIDATES:
        raise BaselineScreenError("too few converged SGD elastic-net candidates")
    selected = min(
        valid_candidates,
        key=lambda item: (
            float(item["mean_inner_log_loss"]),
            float(item["alpha"]),
            float(item["l1_ratio"]),
        ),
    )
    return {
        "alpha": float(selected["alpha"]),
        "l1_ratio": float(selected["l1_ratio"]),
    }, candidate_audit


def fit_ranked_outer_sgd_elastic_net(
    train_x: np.ndarray,
    train_y: np.ndarray,
    sample_weight: np.ndarray,
    candidate_audit: Sequence[Mapping[str, Any]],
    *,
    fold: int,
    seed: int,
) -> tuple[Any, dict[str, Any], dict[str, float], list[dict[str, Any]]]:
    """Fit inner-ranked candidates on outer training data until one converges."""
    ranked = sorted(
        (item for item in candidate_audit if item["valid_all_inner_folds"]),
        key=lambda item: (
            float(item["mean_inner_log_loss"]),
            float(item["alpha"]),
            float(item["l1_ratio"]),
        ),
    )
    if len(ranked) < ELASTIC_NET_MIN_VALID_CANDIDATES:
        raise BaselineScreenError("too few inner-ranked SGD elastic-net candidates")
    attempts: list[dict[str, Any]] = []
    for rank, candidate in enumerate(ranked, start=1):
        selected = {
            "alpha": float(candidate["alpha"]),
            "l1_ratio": float(candidate["l1_ratio"]),
        }
        classifier, state = fit_sgd_elastic_net(
            train_x,
            train_y,
            sample_weight,
            alpha=selected["alpha"],
            l1_ratio=selected["l1_ratio"],
            seed=seed + fold,
        )
        attempts.append(
            {
                "inner_validation_rank": rank,
                "mean_inner_log_loss": float(candidate["mean_inner_log_loss"]),
                "selected": selected,
                "fit_state": state,
            }
        )
        if state["converged_before_ceiling"]:
            return classifier, state, selected, attempts
    raise BaselineScreenError("all inner-ranked outer-training SGD fits failed")


def _fit_fold_models(
    *,
    raw_train: Any,
    gene_ids: Sequence[str],
    train_x: np.ndarray,
    train_y: np.ndarray,
    train_donors: np.ndarray,
    train_row_ids: np.ndarray,
    test_x: np.ndarray,
    inner_assignment: Mapping[str, int],
    fold: int,
    seed: int,
    output: Path,
) -> dict[str, np.ndarray]:
    from sklearn.linear_model import LogisticRegression
    from sklearn.svm import LinearSVC

    weights = np.asarray(
        donor_class_weights(list(train_donors), [ROSTER[index] for index in train_y]),
        dtype=np.float64,
    )
    if set(train_y) != set(range(len(ROSTER))) or set(inner_assignment) != set(train_donors):
        raise BaselineScreenError("training labels or inner donor assignment differ")
    results: dict[str, np.ndarray] = {}
    classifier_train_x, classifier_test_x, classifier_mean, classifier_scale = (
        standardize_projection(train_x, test_x)
    )
    _save_array(output / "classifier_projection_mean.npy", classifier_mean)
    _save_array(output / "classifier_projection_scale.npy", classifier_scale)

    centroid, centroid_state = nearest_centroid_probabilities(
        train_x, train_y, weights, test_x, classes=len(ROSTER)
    )
    results["hvg_pca_nearest_centroid"] = _validate_probabilities(centroid, len(test_x))
    _save_array(output / "nearest_centroid_centroids.npy", centroid_state["centroids"])
    write_json_exclusive(
        output / "nearest_centroid.json",
        {"temperature": centroid_state["temperature"]},
    )

    order = np.argsort(train_row_ids, kind="stable")
    knn = weighted_knn_probabilities(
        train_x[order],
        train_y[order],
        weights[order],
        test_x,
        classes=len(ROSTER),
        neighbors=min(KNN_NEIGHBORS, len(train_x)),
        n_jobs=max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))),
    )
    results["hvg_pca_knn"] = _validate_probabilities(knn, len(test_x))

    logistic = LogisticRegression(
        C=1.0,
        penalty="l2",
        solver="lbfgs",
        tol=1.0e-4,
        max_iter=LOGISTIC_MAX_ITER,
        fit_intercept=True,
        random_state=seed + fold,
    )
    logistic.fit(classifier_train_x, train_y, sample_weight=weights)
    if tuple(logistic.classes_) != tuple(range(len(ROSTER))) or np.any(
        logistic.n_iter_ >= LOGISTIC_MAX_ITER
    ):
        raise BaselineScreenError("hvg_pca_logistic class order or convergence differs")
    results["hvg_pca_logistic"] = _validate_probabilities(
        logistic.predict_proba(classifier_test_x), len(test_x)
    )
    _save_array(output / "hvg_pca_logistic_coef.npy", logistic.coef_)
    _save_array(output / "hvg_pca_logistic_intercept.npy", logistic.intercept_)

    prepared = prepare_inner_projections(
        raw_train=raw_train,
        gene_ids=gene_ids,
        train_y=train_y,
        train_donors=train_donors,
        inner_assignment=inner_assignment,
        fold=fold,
        seed=seed,
    )
    inner_selected_elastic, candidate_audit = select_sgd_elastic_net(
        prepared, train_y, fold=fold, seed=seed
    )
    elastic, elastic_state, selected_elastic, outer_fit_attempts = (
        fit_ranked_outer_sgd_elastic_net(
            classifier_train_x,
            train_y,
            weights,
            candidate_audit,
            fold=fold,
            seed=seed,
        )
    )
    results["hvg_pca_elastic_net"] = _validate_probabilities(
        elastic.predict_proba(classifier_test_x), len(test_x)
    )
    _save_array(output / "hvg_pca_elastic_net_coef.npy", elastic.coef_)
    _save_array(output / "hvg_pca_elastic_net_intercept.npy", elastic.intercept_)
    write_json_exclusive(
        output / "hvg_pca_elastic_net_candidate_audit.json",
        {
            "selection_partition": "donor_grouped_inner_training_only",
            "selection_metric": "donor_class_balanced_multiclass_log_loss",
            "inner_rank_one": inner_selected_elastic,
            "selected": selected_elastic,
            "outer_training_convergence_fallback": outer_fit_attempts,
            "candidates": candidate_audit,
            "test_partition_predicted_during_selection": False,
        },
    )
    logistic_states: dict[str, Any] = {
        "hvg_pca_logistic": {
            "iterations": [int(value) for value in logistic.n_iter_],
            "penalty": "l2",
            "solver": "lbfgs",
            "l1_ratio": None,
            "tolerance": 1.0e-4,
            "max_iter": LOGISTIC_MAX_ITER,
        },
        "hvg_pca_elastic_net": {
            **elastic_state,
            "selection_partition": "donor_grouped_inner_training_only",
            "selection_metric": "donor_class_balanced_multiclass_log_loss",
            "inner_rank_one": inner_selected_elastic,
            "selected": selected_elastic,
            "selected_inner_validation_rank": len(outer_fit_attempts),
            "outer_training_convergence_fallback": outer_fit_attempts,
            "outer_test_partition_used_for_fallback": False,
            "valid_candidates": sum(
                bool(item["valid_all_inner_folds"]) for item in candidate_audit
            ),
        },
    }
    write_json_exclusive(output / "logistic_models.json", logistic_states)

    oof_scores = np.full((len(train_x), len(ROSTER)), np.nan, dtype=np.float64)
    oof_seen = np.zeros(len(train_x), dtype=np.int8)
    for item in prepared:
        validation = np.asarray(item["validation"], dtype=bool)
        fitting = np.asarray(item["fitting"], dtype=bool)
        inner = LinearSVC(
            C=1.0,
            dual="auto",
            max_iter=10_000,
            random_state=seed + fold * 11 + int(item["inner_fold"]),
        )
        inner.fit(
            np.asarray(item["fitting_x"]),
            train_y[fitting],
            sample_weight=np.asarray(item["fitting_weights"]),
        )
        if tuple(inner.classes_) != tuple(range(len(ROSTER))) or int(inner.n_iter_) >= 10_000:
            raise BaselineScreenError("inner linear SVM class order or convergence differs")
        oof_scores[validation] = inner.decision_function(
            np.asarray(item["validation_x"])
        )
        oof_seen[validation] += 1
    if not np.all(np.isfinite(oof_scores)) or not np.all(oof_seen == 1):
        raise BaselineScreenError("SVM cross-fitted calibration is incomplete")
    calibrator = LogisticRegression(
        C=1.0,
        penalty="l2",
        solver="lbfgs",
        max_iter=2_000,
        random_state=seed + fold,
    )
    calibrator.fit(oof_scores, train_y, sample_weight=weights)
    svm = LinearSVC(
        C=1.0,
        dual="auto",
        max_iter=10_000,
        random_state=seed + fold,
    )
    svm.fit(classifier_train_x, train_y, sample_weight=weights)
    if (
        tuple(svm.classes_) != tuple(range(len(ROSTER)))
        or tuple(calibrator.classes_) != tuple(range(len(ROSTER)))
        or int(svm.n_iter_) >= 10_000
        or np.any(calibrator.n_iter_ >= 2_000)
    ):
        raise BaselineScreenError("final SVM or calibrator differs")
    scores = svm.decision_function(classifier_test_x)
    results["hvg_pca_linear_svm"] = _validate_probabilities(
        calibrator.predict_proba(scores), len(test_x)
    )
    for name, value in {
        "svm_coef.npy": svm.coef_,
        "svm_intercept.npy": svm.intercept_,
        "svm_calibration_coef.npy": calibrator.coef_,
        "svm_calibration_intercept.npy": calibrator.intercept_,
    }.items():
        _save_array(output / name, value)
    return results


def run(
    *,
    source: Path,
    split: Path,
    output: Path,
    expected_source_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
) -> None:
    import anndata
    from sklearn.decomposition import PCA

    if output.exists():
        raise BaselineScreenError("refusing to overwrite baseline output")
    if sha256_file(source / "ARTIFACTS.json") != expected_source_artifacts_sha256:
        raise BaselineScreenError("source ARTIFACTS SHA-256 differs")
    if sha256_file(split / "ARTIFACTS.json") != expected_split_artifacts_sha256:
        raise BaselineScreenError("split ARTIFACTS SHA-256 differs")
    source_manifest = json.loads((source / "ARTIFACTS.json").read_text(encoding="utf-8"))
    split_manifest = json.loads((split / "ARTIFACTS.json").read_text(encoding="utf-8"))
    if (
        source_manifest["metadata"].get("subset_id") != DATASET_VIEW_ID
        or source_manifest["metadata"].get("row_count") != 50_000
        or split_manifest["metadata"].get("split_id") != SPLIT_ID
        or split_manifest["metadata"].get("target_labels_used_for_assignment") is not False
        or split_manifest["metadata"].get("sealed_outcomes_read") is not False
    ):
        raise BaselineScreenError("source or split metadata differs")

    selection = _read_tsv(source / "selection.tsv")
    split_rows = _read_tsv(split / "row_outer_folds.tsv")
    inner_rows = _read_tsv(split / "inner_donor_folds.tsv")
    adata = anndata.read_h5ad(source / "resource_atlas_frozen_screen_50000.h5ad")
    row_ids = np.asarray(list(map(str, adata.obs["row_id"])), dtype=str)
    donors = np.asarray(list(map(str, adata.obs["donor_id"])), dtype=str)
    datasets = np.asarray(list(map(str, adata.obs["dataset"])), dtype=str)
    labels = np.asarray(list(map(str, adata.obs["broad_label"])), dtype=str)
    expected_row_ids = [row["row_id"] for row in selection]
    if (
        adata.shape[0] != 50_000
        or row_ids.tolist() != expected_row_ids
        or [row["row_id"] for row in split_rows] != expected_row_ids
        or donors.tolist() != [row["donor_id"] for row in split_rows]
        or datasets.tolist() != [row["dataset"] for row in split_rows]
        or set(labels) != set(ROSTER)
    ):
        raise BaselineScreenError("Atlas, selection, and split row contracts differ")
    outer = np.asarray([int(row["outer_fold"]) for row in split_rows], dtype=np.int8)
    if (
        set(outer) != set(range(5))
        or any(len(set(outer[donors == donor])) != 1 for donor in set(donors))
        or any(len(set(outer[datasets == study])) != 1 for study in set(datasets))
    ):
        raise BaselineScreenError("one donor or study crosses outer folds")
    validate_integer_like_counts(adata.X)

    inner_by_outer: dict[int, dict[str, int]] = defaultdict(dict)
    for row in inner_rows:
        held = int(row["held_outer_fold"])
        donor = row["donor_id"]
        inner_fold = int(row["inner_fold"])
        if held not in range(5) or inner_fold not in range(5):
            raise BaselineScreenError("inner or held fold is outside the frozen roster")
        if donor in inner_by_outer[held]:
            raise BaselineScreenError("inner donor assignment repeats")
        inner_by_outer[held][donor] = inner_fold

    output.mkdir(mode=0o750)
    (output / "folds").mkdir()
    probabilities_by_seed = {
        seed: {
            model_id: np.full((len(row_ids), len(ROSTER)), np.nan, dtype=np.float64)
            for model_id in MODEL_IDS
        }
        for seed in SEEDS
    }
    label_to_index = {label: index for index, label in enumerate(ROSTER)}
    targets = np.asarray([label_to_index[label] for label in labels], dtype=np.int64)
    gene_ids = list(map(str, adata.var["ensembl_id"]))
    fold_receipts: list[dict[str, Any]] = []

    for seed in SEEDS:
        seed_root = output / "folds" / f"seed{seed}"
        seed_root.mkdir()
        for fold in range(5):
            fold_root = seed_root / f"fold{fold}"
            fold_root.mkdir()
            train = np.flatnonzero(outer != fold)
            test = np.flatnonzero(outer == fold)
            if (
                set(inner_by_outer[fold]) != set(donors[train])
                or set(inner_by_outer[fold].values()) != set(range(5))
                or set(datasets[train]) & set(datasets[test])
            ):
                raise BaselineScreenError(
                    "frozen inner assignment or study-held outer firewall differs"
                )
            normalized_train = _log_normalize(adata.X[train], TARGET_SUM)
            selected, feature_records = _select_hvgs(
                normalized_train, gene_ids, n_top=N_TOP_HVG, n_bins=20
            )
            train_dense = normalized_train[:, selected].toarray().astype(np.float32, copy=False)
            normalized_test = _log_normalize(adata.X[test], TARGET_SUM)
            test_dense = normalized_test[:, selected].toarray().astype(np.float32, copy=False)
            pca = PCA(
                n_components=N_PCA,
                svd_solver="randomized",
                whiten=False,
                random_state=seed + fold,
            )
            train_x = pca.fit_transform(train_dense).astype(np.float32, copy=False)
            test_x = pca.transform(test_dense).astype(np.float32, copy=False)
            if not np.all(np.isfinite(train_x)) or not np.all(np.isfinite(test_x)):
                raise BaselineScreenError("PCA embeddings are nonfinite")
            _write_tsv(
                fold_root / "selected_features.tsv",
                ("feature_index", "ensembl_id", "training_mean", "training_variance", "normalized_dispersion"),
                feature_records,
            )
            _save_array(fold_root / "pca_components.npy", pca.components_)
            _save_array(fold_root / "pca_mean.npy", pca.mean_)
            _save_array(fold_root / "pca_explained_variance.npy", pca.explained_variance_)
            fold_predictions = _fit_fold_models(
                raw_train=adata.X[train],
                gene_ids=gene_ids,
                train_x=train_x,
                train_y=targets[train],
                train_donors=donors[train],
                train_row_ids=row_ids[train],
                test_x=test_x,
                inner_assignment=inner_by_outer[fold],
                fold=fold,
                seed=seed,
                output=fold_root,
            )
            for model_id, values in fold_predictions.items():
                probabilities_by_seed[seed][model_id][test] = values
            receipt = {
                "seed": seed,
                "outer_fold": fold,
                "training_rows": len(train),
                "training_donors": len(set(donors[train])),
                "training_studies": sorted(set(datasets[train])),
                "test_rows": len(test),
                "test_donors": len(set(donors[test])),
                "test_studies": sorted(set(datasets[test])),
                "selected_feature_ids_sha256": canonical_sha256(
                    [gene_ids[index] for index in selected]
                ),
                "pca_explained_variance_sum": float(pca.explained_variance_ratio_.sum()),
                "preprocessing_fit_on_outer_training_only": True,
                "classifier_projection_standardization_fit_on_outer_training_only": True,
                "inner_svm_preprocessing_fit_on_inner_training_only": True,
                "inner_svm_projection_standardization_fit_on_inner_training_only": True,
                "inner_svm_weights_fit_on_inner_training_only": True,
                "inner_calibration_group": "donor",
            }
            write_json_exclusive(fold_root / "fold_receipt.json", receipt)
            fold_receipts.append(receipt)

    prediction_root = output / "predictions"
    prediction_root.mkdir()
    fields = (
        "row_id",
        "donor_id",
        "dataset",
        "outer_fold",
        "predicted_class",
        *(f"probability::{label}" for label in ROSTER),
    )
    def write_predictions(path: Path, values: np.ndarray) -> None:
        values = _validate_probabilities(values, len(row_ids))
        predicted = [ROSTER[int(np.argmax(row))] for row in values]
        _write_tsv(
            path,
            fields,
            (
                {
                    "row_id": row_ids[index],
                    "donor_id": donors[index],
                    "dataset": datasets[index],
                    "outer_fold": int(outer[index]),
                    "predicted_class": predicted[index],
                    **{
                        f"probability::{label}": format(float(values[index, class_index]), ".17g")
                        for class_index, label in enumerate(ROSTER)
                    },
                }
                for index in range(len(row_ids))
            ),
        )

    for model_id in MODEL_IDS:
        seed_values = []
        for seed in SEEDS:
            values = _validate_probabilities(
                probabilities_by_seed[seed][model_id], len(row_ids)
            )
            seed_values.append(values)
            write_predictions(
                prediction_root / f"{model_id}--seed-{seed}.tsv", values
            )
        write_predictions(
            prediction_root / f"{model_id}.tsv",
            np.mean(np.stack(seed_values, axis=0), axis=0),
        )

    receipt = {
        "schema_version": "masld-bench-cell-baseline-study-screen-v1",
        "status": "pass_development_predictions",
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": len(row_ids),
        "donors": len(set(donors)),
        "studies": len(set(datasets)),
        "outer_folds": 5,
        "inner_folds": 5,
        "models": list(MODEL_IDS),
        "parameters": {
            "normalization": "library_size_10000_log1p",
            "n_top_hvg": N_TOP_HVG,
            "hvg_mean_bins": 20,
            "pca_components": N_PCA,
            "pca_solver": "randomized",
            "classifier_projection_standardization": "training_only_unweighted_zscore",
            "knn_neighbors": KNN_NEIGHBORS,
            "regularization_c": 1.0,
            "elastic_net_solver": "SGDClassifier_log_loss_elasticnet",
            "elastic_net_alphas": list(ELASTIC_NET_ALPHAS),
            "elastic_net_l1_ratios": list(ELASTIC_NET_L1_RATIOS),
            "elastic_net_tolerance": ELASTIC_NET_TOL,
            "elastic_net_max_iter": ELASTIC_NET_MAX_ITER,
            "elastic_net_n_iter_no_change": ELASTIC_NET_N_ITER_NO_CHANGE,
            "elastic_net_average": True,
            "elastic_net_selection": "donor_grouped_inner_training_only_log_loss",
            "elastic_net_outer_fit_rule": "inner_loss_rank_until_outer_training_convergence",
            "classifier_projection_dtype": "float64",
            "seeds": list(SEEDS),
            "final_prediction": "unweighted_mean_probability_ensemble",
        },
        "folds": fold_receipts,
        "metrics_calculated": False,
        "development_labels_read": ["broad_label"],
        "prediction_tables_contain_observed_labels": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(output / "prediction_receipt.json", receipt)


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True, type=Path)
    value.add_argument("--split", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--expected-source-artifacts-sha256", required=True)
    value.add_argument("--expected-split-artifacts-sha256", required=True)
    return value


def main() -> int:
    arguments = parser().parse_args()
    run(
        source=arguments.source,
        split=arguments.split,
        output=arguments.output,
        expected_source_artifacts_sha256=arguments.expected_source_artifacts_sha256,
        expected_split_artifacts_sha256=arguments.expected_split_artifacts_sha256,
    )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
