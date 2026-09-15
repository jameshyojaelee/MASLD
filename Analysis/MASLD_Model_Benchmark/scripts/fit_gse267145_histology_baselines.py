#!/usr/bin/env python3
"""Fit the seven participant-held GSE267145 task-native baseline families."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from hashlib import sha256
import json
import math
import os
from pathlib import Path
from statistics import mean, stdev
from typing import Any, Callable, Mapping, Sequence
import warnings

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
ELASTIC_LOGISTIC_MAX_ITER = 500_000
FIT_COUNTS: Counter[str] = Counter()
ELASTIC_LOGISTIC_ITERATIONS: list[int] = []


class HistologyBaselineError(RuntimeError):
    """Raised when fitting would violate a lane, fold, or selection requirement."""


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
            raise HistologyBaselineError(f"TSV has no header: {path}")
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


def write_json_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True) + "\n")


def macro_f1(observed: np.ndarray, predicted: np.ndarray, classes: Sequence[int]) -> float:
    scores: list[float] = []
    for class_id in classes:
        tp = int(np.sum((observed == class_id) & (predicted == class_id)))
        fp = int(np.sum((observed != class_id) & (predicted == class_id)))
        fn = int(np.sum((observed == class_id) & (predicted != class_id)))
        denominator = 2 * tp + fp + fn
        scores.append(0.0 if denominator == 0 else 2.0 * tp / denominator)
    return float(np.mean(scores))


def _rank(values: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    return np.asarray(rankdata(values, method="average"), dtype=np.float64)


def spearman(observed: np.ndarray, predicted: np.ndarray) -> float:
    left, right = _rank(observed), _rank(predicted)
    if float(left.std()) == 0.0 or float(right.std()) == 0.0:
        return -1.0
    return float(np.corrcoef(left, right)[0, 1])


def one_standard_error_select(
    candidates: Sequence[Mapping[str, Any]], *, maximize: bool
) -> dict[str, Any]:
    if not candidates:
        raise HistologyBaselineError("one-standard-error selection has no candidates")
    checked: list[dict[str, Any]] = []
    for candidate in candidates:
        scores = [float(value) for value in candidate["fold_scores"]]
        if len(scores) != 4 or not np.all(np.isfinite(scores)):
            continue
        record = dict(candidate)
        record["mean_score"] = mean(scores)
        record["standard_error"] = stdev(scores) / math.sqrt(len(scores))
        checked.append(record)
    if not checked:
        raise HistologyBaselineError("every hyperparameter candidate failed")
    best = min(
        checked,
        key=lambda row: (
            -row["mean_score"] if maximize else row["mean_score"], row["candidate_id"]
        ),
    )
    threshold = (
        best["mean_score"] - best["standard_error"]
        if maximize
        else best["mean_score"] + best["standard_error"]
    )
    if maximize:
        eligible = [row for row in checked if row["mean_score"] >= threshold]
    else:
        eligible = [row for row in checked if row["mean_score"] <= threshold]
    selected = min(eligible, key=lambda row: (tuple(row["complexity"]), row["candidate_id"]))
    selected["one_standard_error_threshold"] = threshold
    selected["one_standard_error_candidate_count"] = len(eligible)
    selected["best_mean_score"] = best["mean_score"]
    selected["best_standard_error"] = best["standard_error"]
    return selected


def select_with_failure_audit(
    candidates: Sequence[Mapping[str, Any]], *, maximize: bool
) -> dict[str, Any]:
    selected = one_standard_error_select(candidates, maximize=maximize)
    failures: Counter[str] = Counter()
    valid_count = 0
    for candidate in candidates:
        scores = candidate.get("fold_scores", [])
        if len(scores) == 4 and np.all(np.isfinite(np.asarray(scores, dtype=float))):
            valid_count += 1
        else:
            failures[str(candidate.get("failure_reason") or "invalid_fold_scores")] += 1
    selected["candidate_failure_count"] = sum(failures.values())
    selected["candidate_failure_reasons"] = dict(sorted(failures.items()))
    selected["candidate_count"] = len(candidates)
    selected["valid_candidate_count"] = valid_count
    return selected


def _log_library_scale(matrix: np.ndarray) -> np.ndarray:
    value = np.asarray(matrix, dtype=np.float64)
    if value.ndim != 2 or not np.all(np.isfinite(value)) or np.any(value < 0):
        raise HistologyBaselineError("molecular matrix is not finite nonnegative")
    totals = value.sum(axis=1)
    if np.any(totals <= 0):
        raise HistologyBaselineError("participant molecular library is empty")
    return np.log1p(value * (1_000_000.0 / totals[:, None])).astype(np.float32)


def _feature_ids(path: Path, field: str) -> list[str]:
    fields, rows = read_tsv(path)
    if field not in fields:
        raise HistologyBaselineError("feature axis field differs")
    values = [row[field] for row in rows]
    if len(values) != len(set(values)) or any(not value for value in values):
        raise HistologyBaselineError("feature identifiers are empty or duplicated")
    return values


class RepresentationCache:
    def __init__(
        self,
        matrix: np.ndarray,
        feature_ids: Sequence[str],
        feature_requests: Sequence[int | str],
        pca_requests: Sequence[int],
        seed: int,
        export_selected_feature_ids: bool = False,
    ) -> None:
        self.matrix = matrix
        self.feature_ids = np.asarray(feature_ids, dtype=str)
        self.feature_requests = tuple(feature_requests)
        self.pca_requests = tuple(int(value) for value in pca_requests)
        self.seed = seed
        self.export_selected_feature_ids = export_selected_feature_ids
        self.cache: dict[tuple[str, str, str, str], dict[str, Any]] = {}

    def get(
        self,
        split_id: str,
        fitting: np.ndarray,
        evaluation: np.ndarray,
        feature_request: int | str,
    ) -> dict[str, Any]:
        from sklearn.decomposition import PCA

        request_id = str(feature_request)
        fitting_hash = sha256(np.asarray(fitting, dtype=np.int64).tobytes()).hexdigest()
        evaluation_hash = sha256(
            np.asarray(evaluation, dtype=np.int64).tobytes()
        ).hexdigest()
        key = (split_id, request_id, fitting_hash, evaluation_hash)
        if key in self.cache:
            return self.cache[key]
        fitting_matrix = self.matrix[fitting]
        minimum_positive = max(3, math.ceil(0.10 * len(fitting)))
        eligible = np.flatnonzero(np.sum(fitting_matrix > 0, axis=0) >= minimum_positive)
        if len(eligible) < 2:
            raise HistologyBaselineError("training partition has fewer than two eligible features")
        variances = np.var(fitting_matrix[:, eligible], axis=0, ddof=0)
        order = np.lexsort((self.feature_ids[eligible], -variances))
        requested = len(eligible) if feature_request == "all_eligible" else int(feature_request)
        count = min(requested, len(eligible))
        selected = eligible[order[:count]]
        fit_dense = np.asarray(fitting_matrix[:, selected], dtype=np.float64)
        eval_dense = np.asarray(self.matrix[evaluation][:, selected], dtype=np.float64)
        center = fit_dense.mean(axis=0)
        scale = fit_dense.std(axis=0)
        scale[scale <= np.finfo(np.float64).eps] = 1.0
        fit_dense = (fit_dense - center) / scale
        eval_dense = (eval_dense - center) / scale
        max_components = min(max(self.pca_requests), len(fitting) - 1, count)
        if max_components < 2:
            raise HistologyBaselineError("training partition cannot support PCA")
        pca_seed = self.seed + int.from_bytes(
            sha256("\0".join(key).encode()).digest()[:4], "big"
        )
        pca = PCA(
            n_components=max_components,
            svd_solver="randomized",
            whiten=False,
            random_state=pca_seed,
        )
        fit_projection = pca.fit_transform(fit_dense)
        eval_projection = pca.transform(eval_dense)
        FIT_COUNTS["pca_fits"] += 1
        selected_ids = [str(value) for value in self.feature_ids[selected]]
        record = {
            "fit": np.asarray(fit_projection, dtype=np.float64),
            "evaluation": np.asarray(eval_projection, dtype=np.float64),
            "actual_feature_count": count,
            "selected_feature_ids_sha256": sha256(
                "\n".join(selected_ids).encode("utf-8")
            ).hexdigest(),
            "max_components": max_components,
            "fitting_index_sha256": fitting_hash,
            "evaluation_index_sha256": evaluation_hash,
        }
        if self.export_selected_feature_ids:
            record["selected_feature_ids"] = selected_ids
        self.cache[key] = record
        return record


def _class_weights(labels: np.ndarray) -> dict[int, float]:
    counts = Counter(int(value) for value in labels)
    return {key: len(labels) / (len(counts) * value) for key, value in counts.items()}


def _logistic(
    x: np.ndarray,
    y: np.ndarray,
    *,
    c_value: float,
    l1_ratio: float,
    seed: int,
):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model = LogisticRegression(
            C=c_value,
            penalty="elasticnet",
            solver="saga",
            l1_ratio=l1_ratio,
            max_iter=ELASTIC_LOGISTIC_MAX_ITER,
            tol=1e-4,
            fit_intercept=True,
            class_weight=_class_weights(y),
            random_state=seed,
        )
        model.fit(x, y)
    FIT_COUNTS["elastic_logistic_fits"] += 1
    ELASTIC_LOGISTIC_ITERATIONS.append(int(np.max(model.n_iter_)))
    if any(isinstance(item.message, ConvergenceWarning) for item in caught):
        raise HistologyBaselineError("elastic-net logistic candidate did not converge")
    return model


def _linear_svm(x: np.ndarray, y: np.ndarray, *, c_value: float, seed: int):
    from sklearn.svm import LinearSVC

    model = LinearSVC(
        C=c_value,
        dual="auto",
        class_weight=_class_weights(y),
        max_iter=20_000,
        random_state=seed,
    )
    model.fit(x, y)
    FIT_COUNTS["linear_svm_fits"] += 1
    if int(model.n_iter_) >= 20_000:
        raise HistologyBaselineError("linear-SVM candidate did not converge")
    return model


def _candidate_id(parameters: Mapping[str, Any]) -> str:
    return sha256(canonical_json(parameters).encode("utf-8")).hexdigest()


def elastic_logistic_solver_diagnostics() -> dict[str, Any]:
    iterations = np.asarray(ELASTIC_LOGISTIC_ITERATIONS, dtype=np.int64)
    return {
        "solver": "saga",
        "tolerance": 1e-4,
        "previous_max_iter": 20_000,
        "intermediate_diagnostic_max_iter": 100_000,
        "max_iter": ELASTIC_LOGISTIC_MAX_ITER,
        "fitted_candidates": int(iterations.size),
        "maximum_iteration_observed": int(iterations.max()) if iterations.size else None,
        "fits_at_iteration_ceiling": int(
            np.count_nonzero(iterations >= ELASTIC_LOGISTIC_MAX_ITER)
        ),
    }


def _temperature_scale(probabilities: np.ndarray, temperature: float) -> np.ndarray:
    if temperature <= 0:
        raise HistologyBaselineError("probability temperature must be positive")
    clipped = np.clip(np.asarray(probabilities, dtype=np.float64), 1e-12, 1.0)
    logits = np.log(clipped) / temperature
    logits -= logits.max(axis=1, keepdims=True)
    scaled = np.exp(logits)
    return scaled / scaled.sum(axis=1, keepdims=True)


def _distance_probabilities(
    kind: str,
    x_fit: np.ndarray,
    y_fit: np.ndarray,
    x_evaluation: np.ndarray,
    *,
    neighbor_count: int,
    temperature: float,
) -> np.ndarray:
    if kind == "nearest_centroid":
        centroids = np.vstack([x_fit[y_fit == label].mean(axis=0) for label in range(3)])
        squared_distance = np.sum(
            (x_evaluation[:, None, :] - centroids[None, :, :]) ** 2, axis=2
        )
        logits = -squared_distance / temperature
        logits -= logits.max(axis=1, keepdims=True)
        probabilities = np.exp(logits)
        FIT_COUNTS["nearest_centroid_fits"] += 1
        return probabilities / probabilities.sum(axis=1, keepdims=True)
    if kind == "knn":
        from sklearn.neighbors import KNeighborsClassifier

        model = KNeighborsClassifier(
            n_neighbors=min(neighbor_count, len(y_fit)), weights="distance", p=2
        )
        model.fit(x_fit, y_fit)
        FIT_COUNTS["knn_fits"] += 1
        probabilities = model.predict_proba(x_evaluation)
        if tuple(int(value) for value in model.classes_) != (0, 1, 2):
            raise HistologyBaselineError("kNN training partition lacks a stage class")
        return _temperature_scale(probabilities, temperature)
    raise HistologyBaselineError(f"unsupported distance classifier: {kind}")


def _fit_multiclass_calibrator(raw: np.ndarray, labels: np.ndarray, seed: int):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model = LogisticRegression(
            C=1.0,
            l1_ratio=0.0,
            max_iter=5_000,
            class_weight=_class_weights(labels),
            random_state=seed,
        )
        model.fit(raw, labels)
    FIT_COUNTS["svm_calibrator_fits"] += 1
    if any(isinstance(item.message, ConvergenceWarning) for item in caught):
        raise HistologyBaselineError("SVM probability calibrator did not converge")
    return model


def _candidate_audit(candidates: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate_id": candidate["candidate_id"],
            "parameters": candidate["parameters"],
            "fold_scores": candidate.get("fold_scores", []),
            "complexity": candidate["complexity"],
            "failure_reason": candidate.get("failure_reason"),
            "outer_training_refit_valid": candidate.get("outer_training_refit_valid", False),
        }
        for candidate in candidates
    ]


def tune_stage_classifier(
    *,
    kind: str,
    cache: RepresentationCache,
    outer_training_indices: np.ndarray,
    inner_assignment: np.ndarray,
    labels: np.ndarray,
    feature_requests: Sequence[int | str],
    pca_requests: Sequence[int],
    c_values: Sequence[float],
    l1_ratios: Sequence[float],
    outer_test_indices: np.ndarray,
    seed: int,
    neighbor_counts: Sequence[int] = (),
    temperatures: Sequence[float] = (1.0,),
    audit_id: str = "stage3",
    audit_callback: Callable[[str, Mapping[str, Any]], None] | None = None,
) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    candidate_outputs: dict[str, dict[str, Any]] = {}
    if kind in {"elastic_net", "svm"}:
        parameter_grid = [
            (feature, pc, c_value, l1, 0, 1.0)
            for feature in feature_requests
            for pc in pca_requests
            for c_value in c_values
            for l1 in (l1_ratios if kind == "elastic_net" else (0.0,))
        ]
    elif kind == "nearest_centroid":
        parameter_grid = [
            (feature, pc, 0.0, 0.0, 0, temperature)
            for feature in feature_requests
            for pc in pca_requests
            for temperature in temperatures
        ]
    elif kind == "knn":
        parameter_grid = [
            (feature, pc, 0.0, 0.0, neighbor_count, temperature)
            for feature in feature_requests
            for pc in pca_requests
            for neighbor_count in neighbor_counts
            for temperature in temperatures
        ]
    else:
        raise HistologyBaselineError(f"unsupported stage classifier: {kind}")
    for feature, pc, c_value, l1, neighbor_count, temperature in parameter_grid:
        scores: list[float] = []
        nonzero: list[int] = []
        actual_counts: list[int] = []
        failure_reason: str | None = None
        oof_raw = np.full((len(outer_training_indices), 3), np.nan, dtype=np.float64)
        full_probabilities: np.ndarray | None = None
        final_rep: Mapping[str, Any] | None = None
        for inner_fold in range(4):
            fit_local = np.flatnonzero(inner_assignment != inner_fold)
            val_local = np.flatnonzero(inner_assignment == inner_fold)
            if set(int(value) for value in labels[fit_local]) != {0, 1, 2}:
                raise HistologyBaselineError("inner stage training lacks a required class")
            fit_global = outer_training_indices[fit_local]
            val_global = outer_training_indices[val_local]
            try:
                rep = cache.get(
                    f"inner{inner_fold}", fit_global, val_global, feature
                )
                components = min(int(pc), int(rep["max_components"]))
                x_fit = rep["fit"][:, :components]
                x_val = rep["evaluation"][:, :components]
                y_fit = labels[fit_local]
                if kind == "elastic_net":
                    model = _logistic(
                        x_fit,
                        y_fit,
                        c_value=float(c_value),
                        l1_ratio=float(l1),
                        seed=seed + inner_fold,
                    )
                    predicted = model.predict(x_val)
                    oof_raw[val_local] = model.predict_proba(x_val)
                    nonzero.append(int(np.count_nonzero(model.coef_)))
                elif kind == "svm":
                    model = _linear_svm(
                        x_fit, y_fit, c_value=float(c_value), seed=seed + inner_fold
                    )
                    predicted = model.predict(x_val)
                    oof_raw[val_local] = model.decision_function(x_val)
                    nonzero.append(int(np.count_nonzero(model.coef_)))
                else:
                    probabilities = _distance_probabilities(
                        kind,
                        x_fit,
                        y_fit,
                        x_val,
                        neighbor_count=int(neighbor_count),
                        temperature=float(temperature),
                    )
                    predicted = np.argmax(probabilities, axis=1)
                    oof_raw[val_local] = probabilities
                    nonzero.append(0)
                scores.append(macro_f1(labels[val_local], predicted, range(3)))
                actual_counts.append(int(rep["actual_feature_count"]))
            except (HistologyBaselineError, ValueError) as error:
                failure_reason = f"{type(error).__name__}:{error}"
                break
        parameters = {
            "feature_request": feature,
            "pca_components": int(pc),
            "c": float(c_value),
            "l1_ratio": float(l1),
            "neighbor_count": int(neighbor_count),
            "temperature": float(temperature),
        }
        candidate_id = _candidate_id(parameters)
        if failure_reason is None:
            try:
                final_rep = cache.get(
                    "outer_final",
                    outer_training_indices,
                    outer_test_indices,
                    feature,
                )
                components = min(int(pc), int(final_rep["max_components"]))
                x_train = final_rep["fit"][:, :components]
                x_test = final_rep["evaluation"][:, :components]
                if kind == "elastic_net":
                    final_model = _logistic(
                        x_train,
                        labels,
                        c_value=float(c_value),
                        l1_ratio=float(l1),
                        seed=seed + 1000,
                    )
                    full_probabilities = final_model.predict_proba(x_test)
                elif kind == "svm":
                    calibrator = _fit_multiclass_calibrator(oof_raw, labels, seed)
                    final_model = _linear_svm(
                        x_train, labels, c_value=float(c_value), seed=seed + 1000
                    )
                    full_probabilities = calibrator.predict_proba(
                        final_model.decision_function(x_test)
                    )
                    oof_raw = calibrator.predict_proba(oof_raw)
                else:
                    full_probabilities = _distance_probabilities(
                        kind,
                        x_train,
                        labels,
                        x_test,
                        neighbor_count=int(neighbor_count),
                        temperature=float(temperature),
                    )
            except (HistologyBaselineError, ValueError) as error:
                failure_reason = f"{type(error).__name__}:{error}"
        candidate = {
                "candidate_id": candidate_id,
                "parameters": parameters,
                "fold_scores": scores if failure_reason is None else [],
                "failure_reason": failure_reason,
                "outer_training_refit_valid": failure_reason is None,
                "complexity": [
                    max(actual_counts) if actual_counts else 10**9,
                    int(pc),
                    mean(nonzero) if nonzero else 10**9,
                    int(neighbor_count),
                    abs(math.log(float(temperature))),
                    float(c_value),
                    float(l1),
                ],
            }
        candidates.append(candidate)
        if failure_reason is None and final_rep is not None and full_probabilities is not None:
            candidate_outputs[candidate_id] = {
                "oof_probabilities": np.asarray(oof_raw, dtype=np.float64),
                "test_probabilities": np.asarray(full_probabilities, dtype=np.float64),
                "final_rep": final_rep,
            }
    selected = select_with_failure_audit(candidates, maximize=True)
    selected["seed_lineage"] = {
        "candidate_inner": "base_seed + inner_fold",
        "candidate_outer_training_refit": "base_seed + 1000",
        "selected_predictions_reuse_exact_candidate_fits": True,
    }
    if audit_callback is not None:
        audit_callback(audit_id, {"selected": selected, "candidates": _candidate_audit(candidates)})
    output = candidate_outputs[selected["candidate_id"]]
    final_rep = output["final_rep"]
    parameters = selected["parameters"]
    return {
        "test_probabilities": output["test_probabilities"],
        "oof_probabilities": output["oof_probabilities"],
        "selected": selected,
        "representation": {
            "feature_request": parameters["feature_request"],
            "pca_components": parameters["pca_components"],
        },
        "final_selected_feature_ids_sha256": final_rep["selected_feature_ids_sha256"],
        "final_selected_feature_ids": final_rep.get("selected_feature_ids"),
    }


def _fit_regression_model(
    x: np.ndarray,
    y: np.ndarray,
    *,
    alpha: float,
    l1_ratio: float,
    seed: int,
):
    from sklearn.exceptions import ConvergenceWarning

    if l1_ratio == 0.0:
        from sklearn.linear_model import Ridge

        model = Ridge(alpha=alpha, fit_intercept=True, solver="lsqr", tol=1e-7)
    else:
        from sklearn.linear_model import ElasticNet

        model = ElasticNet(
            alpha=alpha,
            l1_ratio=l1_ratio,
            max_iter=20_000,
            tol=1e-5,
            random_state=seed,
        )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(x, y)
    FIT_COUNTS["linear_regression_fits"] += 1
    if any(isinstance(item.message, ConvergenceWarning) for item in caught):
        raise HistologyBaselineError("elastic-net regression candidate did not converge")
    return model


def tune_regression(
    *,
    endpoint_id: str,
    cache: RepresentationCache,
    representation: Mapping[str, Any],
    outer_training_indices: np.ndarray,
    outer_test_indices: np.ndarray,
    inner_assignment: np.ndarray,
    target: np.ndarray,
    alphas: Sequence[float],
    l1_ratios: Sequence[float],
    maximize: bool,
    score: Callable[[np.ndarray, np.ndarray], float],
    clip: tuple[float, float],
    seed: int,
    audit_callback: Callable[[str, Mapping[str, Any]], None] | None = None,
) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    candidate_outputs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for alpha in alphas:
        for l1 in l1_ratios:
            scores: list[float] = []
            nonzero: list[int] = []
            failure_reason: str | None = None
            oof = np.full(len(outer_training_indices), np.nan, dtype=np.float64)
            for inner_fold in range(4):
                fit_local = np.flatnonzero(inner_assignment != inner_fold)
                val_local = np.flatnonzero(inner_assignment == inner_fold)
                if len(set(float(value) for value in target[fit_local])) < 2:
                    raise HistologyBaselineError("inner regression training target is constant")
                try:
                    rep = cache.get(
                        f"inner{inner_fold}",
                        outer_training_indices[fit_local],
                        outer_training_indices[val_local],
                        representation["feature_request"],
                    )
                    components = min(representation["pca_components"], rep["max_components"])
                    model = _fit_regression_model(
                        rep["fit"][:, :components],
                        target[fit_local],
                        alpha=float(alpha),
                        l1_ratio=float(l1),
                        seed=seed + inner_fold,
                    )
                    predicted = np.clip(model.predict(rep["evaluation"][:, :components]), *clip)
                    oof[val_local] = predicted
                    scores.append(score(target[val_local], predicted))
                    nonzero.append(int(np.count_nonzero(model.coef_)))
                except (HistologyBaselineError, ValueError) as error:
                    failure_reason = f"{type(error).__name__}:{error}"
                    break
            parameters = {"alpha": float(alpha), "l1_ratio": float(l1)}
            candidate_id = _candidate_id(parameters)
            test: np.ndarray | None = None
            if failure_reason is None:
                try:
                    final_rep = cache.get(
                        "outer_final",
                        outer_training_indices,
                        outer_test_indices,
                        representation["feature_request"],
                    )
                    components = min(
                        representation["pca_components"], final_rep["max_components"]
                    )
                    final = _fit_regression_model(
                        final_rep["fit"][:, :components],
                        target,
                        alpha=float(alpha),
                        l1_ratio=float(l1),
                        seed=seed + 1000,
                    )
                    test = np.clip(
                        final.predict(final_rep["evaluation"][:, :components]), *clip
                    )
                except (HistologyBaselineError, ValueError) as error:
                    failure_reason = f"{type(error).__name__}:{error}"
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "parameters": parameters,
                    "fold_scores": scores if failure_reason is None else [],
                    "failure_reason": failure_reason,
                    "outer_training_refit_valid": failure_reason is None,
                    "complexity": [
                        representation["pca_components"],
                        mean(nonzero) if nonzero else 10**9,
                        -float(alpha),
                        float(l1),
                    ],
                }
            )
            if failure_reason is None and test is not None:
                candidate_outputs[candidate_id] = (oof, np.asarray(test))
    selected = select_with_failure_audit(candidates, maximize=maximize)
    selected["seed_lineage"] = {
        "candidate_inner": "base_seed + inner_fold",
        "candidate_outer_training_refit": "base_seed + 1000",
        "selected_predictions_reuse_exact_candidate_fits": True,
    }
    if audit_callback is not None:
        audit_callback(endpoint_id, {"selected": selected, "candidates": _candidate_audit(candidates)})
    oof, test = candidate_outputs[selected["candidate_id"]]
    return {"oof": oof, "test": np.asarray(test), "selected": selected}


def _pav_decreasing(values: np.ndarray) -> np.ndarray:
    result = np.empty_like(values, dtype=np.float64)
    for row_index, row in enumerate(values):
        blocks = [[float(value), 1] for value in row]
        index = 0
        while index < len(blocks) - 1:
            if blocks[index][0] < blocks[index + 1][0]:
                weight = blocks[index][1] + blocks[index + 1][1]
                average = (
                    blocks[index][0] * blocks[index][1]
                    + blocks[index + 1][0] * blocks[index + 1][1]
                ) / weight
                blocks[index : index + 2] = [[average, weight]]
                index = max(0, index - 1)
            else:
                index += 1
        expanded: list[float] = []
        for value, weight in blocks:
            expanded.extend([value] * weight)
        result[row_index] = np.clip(expanded, 0.0, 1.0)
    return result


def tune_cumulative_fibrosis(
    *,
    cache: RepresentationCache,
    representation: Mapping[str, Any],
    outer_training_indices: np.ndarray,
    outer_test_indices: np.ndarray,
    inner_assignment: np.ndarray,
    fibrosis: np.ndarray,
    c_values: Sequence[float],
    l1_ratios: Sequence[float],
    seed: int,
    audit_callback: Callable[[str, Mapping[str, Any]], None] | None = None,
) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    candidate_outputs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for c_value in c_values:
        for l1 in l1_ratios:
            scores: list[float] = []
            nonzero: list[int] = []
            failure_reason: str | None = None
            oof = np.full(len(outer_training_indices), np.nan, dtype=np.float64)
            for inner_fold in range(4):
                fit_local = np.flatnonzero(inner_assignment != inner_fold)
                val_local = np.flatnonzero(inner_assignment == inner_fold)
                for threshold in (1, 2, 3):
                    binary = (fibrosis[fit_local] >= threshold).astype(np.int64)
                    if set(binary) != {0, 1}:
                        raise HistologyBaselineError("inner cumulative threshold lacks a class")
                try:
                    rep = cache.get(
                        f"inner{inner_fold}",
                        outer_training_indices[fit_local],
                        outer_training_indices[val_local],
                        representation["feature_request"],
                    )
                    components = min(representation["pca_components"], rep["max_components"])
                    probabilities: list[np.ndarray] = []
                    for threshold in (1, 2, 3):
                        binary = (fibrosis[fit_local] >= threshold).astype(np.int64)
                        model = _logistic(
                            rep["fit"][:, :components],
                            binary,
                            c_value=float(c_value),
                            l1_ratio=float(l1),
                            seed=seed + inner_fold * 10 + threshold,
                        )
                        probabilities.append(model.predict_proba(rep["evaluation"][:, :components])[:, 1])
                        nonzero.append(int(np.count_nonzero(model.coef_)))
                    expected = _pav_decreasing(np.column_stack(probabilities)).sum(axis=1)
                    oof[val_local] = expected
                    scores.append(float(np.mean(np.abs(fibrosis[val_local] - expected))))
                except (HistologyBaselineError, ValueError) as error:
                    failure_reason = f"{type(error).__name__}:{error}"
                    break
            parameters = {"c": float(c_value), "l1_ratio": float(l1)}
            candidate_id = _candidate_id(parameters)
            test: np.ndarray | None = None
            if failure_reason is None:
                try:
                    final_rep = cache.get(
                        "outer_final",
                        outer_training_indices,
                        outer_test_indices,
                        representation["feature_request"],
                    )
                    components = min(
                        representation["pca_components"], final_rep["max_components"]
                    )
                    columns: list[np.ndarray] = []
                    for threshold in (1, 2, 3):
                        model = _logistic(
                            final_rep["fit"][:, :components],
                            (fibrosis >= threshold).astype(np.int64),
                            c_value=float(c_value),
                            l1_ratio=float(l1),
                            seed=seed + 1000 + threshold,
                        )
                        columns.append(
                            model.predict_proba(
                                final_rep["evaluation"][:, :components]
                            )[:, 1]
                        )
                    test = _pav_decreasing(np.column_stack(columns)).sum(axis=1)
                except (HistologyBaselineError, ValueError) as error:
                    failure_reason = f"{type(error).__name__}:{error}"
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "parameters": parameters,
                    "fold_scores": scores if failure_reason is None else [],
                    "failure_reason": failure_reason,
                    "outer_training_refit_valid": failure_reason is None,
                    "complexity": [
                        representation["pca_components"],
                        mean(nonzero) if nonzero else 10**9,
                        float(c_value),
                        float(l1),
                    ],
                }
            )
            if failure_reason is None and test is not None:
                candidate_outputs[candidate_id] = (oof, np.asarray(test))
    selected = select_with_failure_audit(candidates, maximize=False)
    selected["seed_lineage"] = {
        "candidate_inner": "base_seed + 10 * inner_fold + threshold",
        "candidate_outer_training_refit": "base_seed + 1000 + threshold",
        "selected_predictions_reuse_exact_candidate_fits": True,
    }
    if audit_callback is not None:
        audit_callback(
            "fibrosis_cumulative",
            {"selected": selected, "candidates": _candidate_audit(candidates)},
        )
    oof, test = candidate_outputs[selected["candidate_id"]]
    return {"oof": oof, "test": test, "selected": selected}


def tune_group_classifier(
    *,
    cache: RepresentationCache,
    representation: Mapping[str, Any],
    outer_training_indices: np.ndarray,
    outer_test_indices: np.ndarray,
    inner_assignment: np.ndarray,
    groups: np.ndarray,
    c_values: Sequence[float],
    l1_ratios: Sequence[float],
    seed: int,
    audit_callback: Callable[[str, Mapping[str, Any]], None] | None = None,
) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    candidate_outputs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for c_value in c_values:
        for l1 in l1_ratios:
            scores: list[float] = []
            nonzero: list[int] = []
            failure_reason: str | None = None
            oof = np.full((len(outer_training_indices), 3), np.nan, dtype=np.float64)
            for inner_fold in range(4):
                fit_local = np.flatnonzero(inner_assignment != inner_fold)
                val_local = np.flatnonzero(inner_assignment == inner_fold)
                if set(int(value) for value in groups[fit_local]) != {0, 1, 2}:
                    raise HistologyBaselineError("inner fibrosis-group training lacks a class")
                try:
                    rep = cache.get(
                        f"inner{inner_fold}",
                        outer_training_indices[fit_local],
                        outer_training_indices[val_local],
                        representation["feature_request"],
                    )
                    components = min(representation["pca_components"], rep["max_components"])
                    model = _logistic(
                        rep["fit"][:, :components],
                        groups[fit_local],
                        c_value=float(c_value),
                        l1_ratio=float(l1),
                        seed=seed + inner_fold,
                    )
                    predicted = model.predict(rep["evaluation"][:, :components])
                    oof[val_local] = model.predict_proba(
                        rep["evaluation"][:, :components]
                    )
                    scores.append(macro_f1(groups[val_local], predicted, range(3)))
                    nonzero.append(int(np.count_nonzero(model.coef_)))
                except (HistologyBaselineError, ValueError) as error:
                    failure_reason = f"{type(error).__name__}:{error}"
                    break
            parameters = {"c": float(c_value), "l1_ratio": float(l1)}
            candidate_id = _candidate_id(parameters)
            test: np.ndarray | None = None
            if failure_reason is None:
                try:
                    final_rep = cache.get(
                        "outer_final",
                        outer_training_indices,
                        outer_test_indices,
                        representation["feature_request"],
                    )
                    components = min(
                        representation["pca_components"], final_rep["max_components"]
                    )
                    final = _logistic(
                        final_rep["fit"][:, :components],
                        groups,
                        c_value=float(c_value),
                        l1_ratio=float(l1),
                        seed=seed + 1000,
                    )
                    test = final.predict_proba(
                        final_rep["evaluation"][:, :components]
                    )
                except (HistologyBaselineError, ValueError) as error:
                    failure_reason = f"{type(error).__name__}:{error}"
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "parameters": parameters,
                    "fold_scores": scores if failure_reason is None else [],
                    "failure_reason": failure_reason,
                    "outer_training_refit_valid": failure_reason is None,
                    "complexity": [
                        representation["pca_components"],
                        mean(nonzero) if nonzero else 10**9,
                        float(c_value),
                        float(l1),
                    ],
                }
            )
            if failure_reason is None and test is not None:
                candidate_outputs[candidate_id] = (oof, np.asarray(test))
    selected = select_with_failure_audit(candidates, maximize=True)
    selected["seed_lineage"] = {
        "candidate_inner": "base_seed + inner_fold",
        "candidate_outer_training_refit": "base_seed + 1000",
        "selected_predictions_reuse_exact_candidate_fits": True,
    }
    if audit_callback is not None:
        audit_callback(
            "fibrosis_group3",
            {"selected": selected, "candidates": _candidate_audit(candidates)},
        )
    oof, test = candidate_outputs[selected["candidate_id"]]
    return {"oof": oof, "test": test, "selected": selected}


def fit_secondary_heads(
    *,
    cache: RepresentationCache,
    representation: Mapping[str, Any],
    outer_training_indices: np.ndarray,
    outer_test_indices: np.ndarray,
    inner_assignment: np.ndarray,
    nas: np.ndarray,
    fibrosis: np.ndarray,
    groups: np.ndarray,
    surface: Mapping[str, Any],
    seed: int,
    audit_callback: Callable[[str, Mapping[str, Any]], None] | None = None,
) -> dict[str, Any]:
    regression = surface["shared_task_native_heads"]["fibrosis_exact_regression"]
    stage_grid = surface["lanes"]["rna_only"]["baselines"]["rna_hvg_pca_elastic_net"]
    nas_result = tune_regression(
        endpoint_id="nash_crn_component_sum",
        cache=cache,
        representation=representation,
        outer_training_indices=outer_training_indices,
        outer_test_indices=outer_test_indices,
        inner_assignment=inner_assignment,
        target=nas,
        alphas=regression["alpha"],
        l1_ratios=regression["l1_ratio"],
        maximize=True,
        score=spearman,
        clip=(0.0, 8.0),
        seed=seed + 10_000,
        audit_callback=audit_callback,
    )
    fibrosis_regression = tune_regression(
        endpoint_id="fibrosis_exact_regression",
        cache=cache,
        representation=representation,
        outer_training_indices=outer_training_indices,
        outer_test_indices=outer_test_indices,
        inner_assignment=inner_assignment,
        target=fibrosis,
        alphas=regression["alpha"],
        l1_ratios=regression["l1_ratio"],
        maximize=False,
        score=lambda observed, predicted: float(np.mean(np.abs(observed - predicted))),
        clip=(0.0, 3.0),
        seed=seed + 20_000,
        audit_callback=audit_callback,
    )
    cumulative = tune_cumulative_fibrosis(
        cache=cache,
        representation=representation,
        outer_training_indices=outer_training_indices,
        outer_test_indices=outer_test_indices,
        inner_assignment=inner_assignment,
        fibrosis=fibrosis,
        c_values=stage_grid["inverse_regularization"],
        l1_ratios=stage_grid["l1_ratio"],
        seed=seed + 30_000,
        audit_callback=audit_callback,
    )
    group = tune_group_classifier(
        cache=cache,
        representation=representation,
        outer_training_indices=outer_training_indices,
        outer_test_indices=outer_test_indices,
        inner_assignment=inner_assignment,
        groups=groups,
        c_values=stage_grid["inverse_regularization"],
        l1_ratios=stage_grid["l1_ratio"],
        seed=seed + 40_000,
        audit_callback=audit_callback,
    )
    return {
        "nas": nas_result,
        "fibrosis_regression": fibrosis_regression,
        "fibrosis_cumulative": cumulative,
        "fibrosis_group": group,
    }


def _prediction_payload(stage: Mapping[str, Any], secondary: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "stage_oof": stage["oof_probabilities"],
        "stage_test": stage["test_probabilities"],
        "nas_oof": secondary["nas"]["oof"],
        "nas_test": secondary["nas"]["test"],
        "fibrosis_cumulative_oof": secondary["fibrosis_cumulative"]["oof"],
        "fibrosis_cumulative_test": secondary["fibrosis_cumulative"]["test"],
        "fibrosis_regression_oof": secondary["fibrosis_regression"]["oof"],
        "fibrosis_regression_test": secondary["fibrosis_regression"]["test"],
        "group_oof": secondary["fibrosis_group"]["oof"],
        "group_test": secondary["fibrosis_group"]["test"],
    }


def _selection_receipt(stage: Mapping[str, Any], secondary: Mapping[str, Any]) -> dict[str, Any]:
    receipt = {
        "stage3": stage["selected"],
        "representation": stage["representation"],
        "final_selected_feature_ids_sha256": stage["final_selected_feature_ids_sha256"],
        "nash_crn_component_sum": secondary["nas"]["selected"],
        "fibrosis_exact_regression": secondary["fibrosis_regression"]["selected"],
        "fibrosis_cumulative": secondary["fibrosis_cumulative"]["selected"],
        "fibrosis_group3": secondary["fibrosis_group"]["selected"],
        "one_standard_error_applied_to_every_selection": True,
        "source_stage5_used_for_selection": False,
        "recorded_sex_used_for_selection": False,
    }
    if stage.get("final_selected_feature_ids") is not None:
        receipt["final_selected_feature_ids"] = stage["final_selected_feature_ids"]
    return receipt


def _candidate_failure_audit(value: Any) -> tuple[int, Counter[str]]:
    total = 0
    reasons: Counter[str] = Counter()
    if isinstance(value, Mapping):
        if "candidate_failure_count" in value:
            total += int(value["candidate_failure_count"])
            reasons.update(
                {
                    str(key): int(count)
                    for key, count in value.get("candidate_failure_reasons", {}).items()
                }
            )
        else:
            for child in value.values():
                child_total, child_reasons = _candidate_failure_audit(child)
                total += child_total
                reasons.update(child_reasons)
    return total, reasons


def _blend_select(
    left: np.ndarray,
    right: np.ndarray,
    target: np.ndarray,
    inner_assignment: np.ndarray,
    weights: Sequence[float],
    *,
    kind: str,
) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    for weight in weights:
        blended = weight * left + (1.0 - weight) * right
        scores: list[float] = []
        for inner_fold in range(4):
            validation = inner_assignment == inner_fold
            if kind == "stage":
                scores.append(macro_f1(target[validation], np.argmax(blended[validation], axis=1), range(3)))
            elif kind == "group":
                scores.append(macro_f1(target[validation], np.argmax(blended[validation], axis=1), range(3)))
            elif kind == "nas":
                scores.append(spearman(target[validation], blended[validation]))
            else:
                scores.append(float(np.mean(np.abs(target[validation] - blended[validation]))))
        parameters = {"left_weight": float(weight)}
        candidates.append(
            {
                "candidate_id": _candidate_id(parameters),
                "parameters": parameters,
                "fold_scores": scores,
                "complexity": [abs(float(weight) - 0.5), canonical_json(parameters)],
            }
        )
    return select_with_failure_audit(
        candidates, maximize=kind in {"stage", "group", "nas"}
    )


def _late_fusion(
    left: Mapping[str, np.ndarray],
    right: Mapping[str, np.ndarray],
    *,
    stage: np.ndarray,
    nas: np.ndarray,
    fibrosis: np.ndarray,
    groups: np.ndarray,
    inner_assignment: np.ndarray,
    weights: Sequence[float],
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    endpoints = {
        "stage": ("stage_oof", "stage_test", stage),
        "nas": ("nas_oof", "nas_test", nas),
        "fibrosis_cumulative": ("fibrosis_cumulative_oof", "fibrosis_cumulative_test", fibrosis),
        "fibrosis_regression": ("fibrosis_regression_oof", "fibrosis_regression_test", fibrosis),
        "group": ("group_oof", "group_test", groups),
    }
    result: dict[str, np.ndarray] = {}
    selections: dict[str, Any] = {}
    for endpoint, (oof_key, test_key, target) in endpoints.items():
        kind = endpoint if endpoint in {"stage", "nas", "group"} else "mae"
        selected = _blend_select(
            left[oof_key], right[oof_key], target, inner_assignment, weights, kind=kind
        )
        weight = selected["parameters"]["left_weight"]
        result[oof_key] = weight * left[oof_key] + (1.0 - weight) * right[oof_key]
        result[test_key] = weight * left[test_key] + (1.0 - weight) * right[test_key]
        selections[endpoint] = selected
    selections["one_standard_error_applied_to_every_selection"] = True
    selections["source_stage5_used_for_selection"] = False
    selections["recorded_sex_used_for_selection"] = False
    return result, selections


class BlockRepresentationCache:
    def __init__(
        self,
        left: RepresentationCache,
        right: RepresentationCache,
        left_representation: Mapping[str, Any],
        right_representation: Mapping[str, Any],
    ) -> None:
        self.left = left
        self.right = right
        self.left_representation = left_representation
        self.right_representation = right_representation
        self.cache: dict[tuple[str, str, str], dict[str, Any]] = {}
        self.feature_requests = ("locked_blocks",)
        self.pca_requests = (
            int(left_representation["pca_components"])
            + int(right_representation["pca_components"]),
        )

    def get(
        self,
        split_id: str,
        fitting: np.ndarray,
        evaluation: np.ndarray,
        feature_request: int | str,
    ) -> dict[str, Any]:
        if feature_request != "locked_blocks":
            raise HistologyBaselineError("block feature request differs")
        fitting_hash = sha256(np.asarray(fitting, dtype=np.int64).tobytes()).hexdigest()
        evaluation_hash = sha256(
            np.asarray(evaluation, dtype=np.int64).tobytes()
        ).hexdigest()
        key = (split_id, fitting_hash, evaluation_hash)
        if key in self.cache:
            return self.cache[key]
        left = self.left.get(
            split_id,
            fitting,
            evaluation,
            self.left_representation["feature_request"],
        )
        right = self.right.get(
            split_id,
            fitting,
            evaluation,
            self.right_representation["feature_request"],
        )
        left_pc = min(self.left_representation["pca_components"], left["max_components"])
        right_pc = min(self.right_representation["pca_components"], right["max_components"])
        left_fit, left_eval = left["fit"][:, :left_pc], left["evaluation"][:, :left_pc]
        right_fit, right_eval = right["fit"][:, :right_pc], right["evaluation"][:, :right_pc]
        left_norm = max(float(np.linalg.norm(left_fit)), np.finfo(float).eps)
        right_norm = max(float(np.linalg.norm(right_fit)), np.finfo(float).eps)
        fit = np.column_stack((left_fit / left_norm, right_fit / right_norm))
        evaluation_values = np.column_stack((left_eval / left_norm, right_eval / right_norm))
        record = {
            "fit": fit,
            "evaluation": evaluation_values,
            "actual_feature_count": left["actual_feature_count"] + right["actual_feature_count"],
            "selected_feature_ids_sha256": sha256(
                (left["selected_feature_ids_sha256"] + right["selected_feature_ids_sha256"]).encode()
            ).hexdigest(),
            "max_components": fit.shape[1],
            "fitting_index_sha256": fitting_hash,
            "evaluation_index_sha256": evaluation_hash,
        }
        self.cache[key] = record
        return record


def _fit_outer_seed(
    *,
    outer_fold: int,
    seed: int,
    rna: np.ndarray,
    h3: np.ndarray,
    rna_ids: Sequence[str],
    h3_ids: Sequence[str],
    participant_ids: Sequence[str],
    outer_values: np.ndarray,
    fit_view: Path,
    surface: Mapping[str, Any],
    output: Path,
) -> dict[str, Any]:
    FIT_COUNTS.clear()
    ELASTIC_LOGISTIC_ITERATIONS.clear()
    fold_root = output / f"outer_{outer_fold}" / f"seed_{seed}"
    fold_root.mkdir(parents=True)
    _, training_rows = read_tsv(fit_view / f"outer_{outer_fold}/training_endpoints.tsv")
    _, query_rows = read_tsv(fit_view / f"outer_{outer_fold}/query_participants.tsv")
    lookup = {participant: index for index, participant in enumerate(participant_ids)}
    training_indices = np.asarray([lookup[row["participant_id"]] for row in training_rows], dtype=int)
    test_indices = np.asarray([lookup[row["participant_id"]] for row in query_rows], dtype=int)
    if np.any(outer_values[training_indices] == outer_fold) or np.any(outer_values[test_indices] != outer_fold):
        raise HistologyBaselineError("fit view crosses the outer fold")
    inner = np.asarray([int(row["inner_validation_fold"]) for row in training_rows], dtype=int)
    stage = np.asarray([STAGE3.index(row["stage3"]) for row in training_rows], dtype=int)
    nas = np.asarray([int(row["nash_crn_component_sum"]) for row in training_rows], dtype=float)
    fibrosis = np.asarray([int(row["fibrosis"]) for row in training_rows], dtype=float)
    groups = np.asarray([FIBROSIS_GROUP3.index(row["fibrosis_group3"]) for row in training_rows], dtype=int)
    rna_lane = surface["lanes"]["rna_only"]
    h3_lane = surface["lanes"]["h3k27ac_only"]
    rna_cache = RepresentationCache(
        rna, rna_ids, rna_lane["feature_screen"]["candidate_counts"],
        rna_lane["baselines"]["rna_hvg_pca_elastic_net"]["pca_components"], seed + outer_fold * 100_000,
        export_selected_feature_ids=True,
    )
    h3_cache = RepresentationCache(
        h3, h3_ids, h3_lane["feature_screen"]["candidate_counts"],
        h3_lane["baselines"]["h3_variance_pca_elastic_net"]["pca_components"], seed + outer_fold * 100_000 + 50_000,
    )
    results: dict[str, dict[str, np.ndarray]] = {}
    receipts: dict[str, Any] = {}
    model_objects: dict[str, tuple[RepresentationCache, dict[str, Any]]] = {}
    for modality, cache, lane, pairs in (
        (
            "rna", rna_cache, rna_lane,
            (
                ("rna_hvg_pca_elastic_net", "elastic_net"),
                ("rna_hvg_pca_linear_svm", "svm"),
                ("rna_hvg_pca_nearest_centroid", "nearest_centroid"),
                ("rna_hvg_pca_knn", "knn"),
            ),
        ),
        (
            "h3", h3_cache, h3_lane,
            (
                ("h3_variance_pca_elastic_net", "elastic_net"),
                ("h3_variance_pca_linear_svm", "svm"),
                ("h3_variance_pca_nearest_centroid", "nearest_centroid"),
                ("h3_variance_pca_knn", "knn"),
            ),
        ),
    ):
        for model_id, kind in pairs:
            spec = lane["baselines"][model_id]
            def audit_callback(
                endpoint_id: str,
                payload: Mapping[str, Any],
                *,
                current_model_id: str = model_id,
            ) -> None:
                write_json_exclusive(
                    fold_root / f"candidate_audit--{current_model_id}--{endpoint_id}.json",
                    {
                        "model_id": current_model_id,
                        "endpoint_id": endpoint_id,
                        **payload,
                    },
                )

            stage_result = tune_stage_classifier(
                kind=kind,
                cache=cache,
                outer_training_indices=training_indices,
                inner_assignment=inner,
                labels=stage,
                feature_requests=lane["feature_screen"]["candidate_counts"],
                pca_requests=spec["pca_components"],
                c_values=spec.get("inverse_regularization", []),
                l1_ratios=spec.get("l1_ratio", [0.0]),
                outer_test_indices=test_indices,
                seed=seed + outer_fold * 1000 + len(receipts) * 100,
                neighbor_counts=spec.get("neighbor_counts", []),
                temperatures=spec.get("probability_temperature", [1.0]),
                audit_id="stage3",
                audit_callback=audit_callback,
            )
            secondary = fit_secondary_heads(
                cache=cache,
                representation=stage_result["representation"],
                outer_training_indices=training_indices,
                outer_test_indices=test_indices,
                inner_assignment=inner,
                nas=nas,
                fibrosis=fibrosis,
                groups=groups,
                surface=surface,
                seed=seed + outer_fold * 10_000 + len(receipts) * 1000,
                audit_callback=audit_callback,
            )
            results[model_id] = _prediction_payload(stage_result, secondary)
            receipts[model_id] = _selection_receipt(stage_result, secondary)
            model_objects[model_id] = (cache, stage_result["representation"])
    block_spec = surface["lanes"]["observed_pair_multimodal"]["baselines"]["block_pca_elastic_net"]
    block_cache = BlockRepresentationCache(
        rna_cache,
        h3_cache,
        model_objects["rna_hvg_pca_elastic_net"][1],
        model_objects["h3_variance_pca_elastic_net"][1],
    )
    block_stage = tune_stage_classifier(
        kind="elastic_net",
        cache=block_cache,  # type: ignore[arg-type]
        outer_training_indices=training_indices,
        inner_assignment=inner,
        labels=stage,
        feature_requests=["locked_blocks"],
        pca_requests=[sum(block_cache.pca_requests)],
        c_values=block_spec["inverse_regularization"],
        l1_ratios=block_spec["l1_ratio"],
        outer_test_indices=test_indices,
        seed=seed + outer_fold * 1000 + 800,
        audit_id="stage3",
        audit_callback=lambda endpoint_id, payload: write_json_exclusive(
            fold_root / f"candidate_audit--block_pca_elastic_net--{endpoint_id}.json",
            {"model_id": "block_pca_elastic_net", "endpoint_id": endpoint_id, **payload},
        ),
    )
    block_secondary = fit_secondary_heads(
        cache=block_cache,  # type: ignore[arg-type]
        representation=block_stage["representation"],
        outer_training_indices=training_indices,
        outer_test_indices=test_indices,
        inner_assignment=inner,
        nas=nas,
        fibrosis=fibrosis,
        groups=groups,
        surface=surface,
        seed=seed + outer_fold * 10_000 + 8000,
        audit_callback=lambda endpoint_id, payload: write_json_exclusive(
            fold_root / f"candidate_audit--block_pca_elastic_net--{endpoint_id}.json",
            {"model_id": "block_pca_elastic_net", "endpoint_id": endpoint_id, **payload},
        ),
    )
    results["block_pca_elastic_net"] = _prediction_payload(block_stage, block_secondary)
    receipts["block_pca_elastic_net"] = _selection_receipt(block_stage, block_secondary)
    fusion_spec = surface["lanes"]["observed_pair_multimodal"]["baselines"]["calibrated_late_fusion"]
    fused, fusion_receipt = _late_fusion(
        results["rna_hvg_pca_elastic_net"],
        results["h3_variance_pca_elastic_net"],
        stage=stage,
        nas=nas,
        fibrosis=fibrosis,
        groups=groups,
        inner_assignment=inner,
        weights=fusion_spec["rna_weight_grid"],
    )
    results["calibrated_late_fusion"] = fused
    receipts["calibrated_late_fusion"] = fusion_receipt
    stage_counts = np.bincount(stage, minlength=3).astype(float)
    stage_probability = stage_counts / stage_counts.sum()
    group_counts = np.bincount(groups, minlength=3).astype(float)
    group_probability = group_counts / group_counts.sum()
    prevalence = {
        "stage_test": np.repeat(stage_probability[None, :], len(test_indices), axis=0),
        "nas_test": np.full(len(test_indices), float(np.median(nas))),
        "fibrosis_cumulative_test": np.full(len(test_indices), float(np.median(fibrosis))),
        "fibrosis_regression_test": np.full(len(test_indices), float(np.median(fibrosis))),
        "group_test": np.repeat(group_probability[None, :], len(test_indices), axis=0),
    }
    results["training_stage_distribution"] = prevalence
    receipts["training_stage_distribution"] = {
        "training_only_class_proportions": True,
        "training_only_continuous_medians": True,
        "one_standard_error_not_applicable_no_hyperparameters": True,
        "source_stage5_used_for_selection": False,
        "recorded_sex_used_for_selection": False,
    }
    candidate_failure_count, candidate_failure_reasons = _candidate_failure_audit(receipts)
    receipts["fit_instrumentation"] = {
        "fit_counts": dict(sorted(FIT_COUNTS.items())),
        "elastic_logistic_solver": elastic_logistic_solver_diagnostics(),
        "cache_identity": "split_id_plus_fitting_and_evaluation_index_sha256",
        "invalid_candidate_count": candidate_failure_count,
        "invalid_candidate_reasons": dict(sorted(candidate_failure_reasons.items())),
    }
    (fold_root / "selection_receipts.json").write_text(
        json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        "outer_fold": outer_fold,
        "seed": seed,
        "test_indices": test_indices,
        "results": {model_id: {key: value for key, value in payload.items() if key.endswith("_test")} for model_id, payload in results.items()},
        "fit_counts": dict(sorted(FIT_COUNTS.items())),
        "elastic_logistic_solver": elastic_logistic_solver_diagnostics(),
        "invalid_candidate_count": candidate_failure_count,
        "invalid_candidate_reasons": dict(sorted(candidate_failure_reasons.items())),
    }


def _fit_outer_bundle(
    *,
    outer_fold: int,
    seeds: Sequence[int],
    molecular: Path,
    folds: Path,
    fit_views: Path,
    surface: Mapping[str, Any],
    output: Path,
) -> list[dict[str, Any]]:
    _, participants = read_tsv(molecular / "participant_axis.tsv")
    _, fold_rows = read_tsv(folds / "participant_outer_folds.tsv")
    participant_ids = [row["participant_id"] for row in participants]
    outer = np.asarray([int(row["outer_fold"]) for row in fold_rows], dtype=int)
    rna_raw = np.load(molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    h3_raw = np.load(molecular / "h3k27ac_counts.npy", mmap_mode="r", allow_pickle=False)
    rna = _log_library_scale(rna_raw)
    h3 = _log_library_scale(h3_raw)
    rna_ids = _feature_ids(molecular / "rna_feature_axis.tsv", "stable_gene_id")
    h3_ids = _feature_ids(
        molecular / "h3k27ac_feature_axis.tsv", "opaque_source_feature_key"
    )
    if rna.shape != (99, 42163) or h3.shape != (99, 96460):
        raise HistologyBaselineError("molecular fixture shapes differ")
    results: list[dict[str, Any]] = []
    for seed in seeds:
        try:
            results.append(
                _fit_outer_seed(
                    outer_fold=outer_fold,
                    seed=seed,
                    rna=rna,
                    h3=h3,
                    rna_ids=rna_ids,
                    h3_ids=h3_ids,
                    participant_ids=participant_ids,
                    outer_values=outer,
                    fit_view=fit_views,
                    surface=surface,
                    output=output,
                )
            )
        except Exception as error:
            write_fit_failure_receipt(output, outer_fold, seed, error)
            raise
    return results


def write_fit_failure_receipt(
    output: Path, outer_fold: int, seed: int, error: Exception
) -> Path:
    fold_root = output / f"outer_{outer_fold}" / f"seed_{seed}"
    fold_root.mkdir(parents=True, exist_ok=True)
    candidate_audits = sorted(fold_root.glob("candidate_audit--*.json"))
    failure: dict[str, Any] = {
        "schema_version": "masld-bench-gse267145-fit-failure-v1",
        "outer_fold": outer_fold,
        "seed": seed,
        "failure_class": type(error).__name__,
        "failure_message": str(error),
        "fit_counts_at_failure": dict(sorted(FIT_COUNTS.items())),
        "elastic_logistic_solver_at_failure": elastic_logistic_solver_diagnostics(),
        "candidate_audit_files_persisted": len(candidate_audits),
        "candidate_audit_file_names": [path.name for path in candidate_audits],
        "prediction_files_written": 0,
        "status": "failed_closed",
    }
    path = fold_root / "failure_receipt.json"
    write_json_exclusive(path, failure)
    return path


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


def _prediction_rows(
    participant_ids: Sequence[str], outer: np.ndarray, payload: Mapping[str, np.ndarray]
) -> list[dict[str, Any]]:
    stage_prob = np.asarray(payload["stage_test"], dtype=float)
    group_prob = np.asarray(payload["group_test"], dtype=float)
    participant_count = len(participant_ids)
    if (
        stage_prob.shape != (participant_count, 3)
        or group_prob.shape != (participant_count, 3)
        or len(outer) != participant_count
        or not np.allclose(stage_prob.sum(axis=1), 1.0, atol=1e-6)
        or not np.allclose(group_prob.sum(axis=1), 1.0, atol=1e-6)
    ):
        raise HistologyBaselineError("assembled prediction probabilities differ")
    rows: list[dict[str, Any]] = []
    for index, participant in enumerate(participant_ids):
        rows.append(
            {
                "participant_id": participant,
                "outer_fold": int(outer[index]),
                "probability_NOR": format(stage_prob[index, 0], ".17g"),
                "probability_NAFL": format(stage_prob[index, 1], ".17g"),
                "probability_NASH": format(stage_prob[index, 2], ".17g"),
                "predicted_stage3": STAGE3[int(np.argmax(stage_prob[index]))],
                "predicted_nash_crn_component_sum": format(float(payload["nas_test"][index]), ".17g"),
                "predicted_fibrosis_cumulative_expected": format(float(payload["fibrosis_cumulative_test"][index]), ".17g"),
                "predicted_fibrosis_regression": format(float(payload["fibrosis_regression_test"][index]), ".17g"),
                "probability_fibrosis_F0": format(group_prob[index, 0], ".17g"),
                "probability_fibrosis_F1": format(group_prob[index, 1], ".17g"),
                "probability_fibrosis_F2_3": format(group_prob[index, 2], ".17g"),
                "predicted_fibrosis_group3": FIBROSIS_GROUP3[int(np.argmax(group_prob[index]))],
            }
        )
    return rows


def fit_preflight(
    *,
    molecular: Path,
    folds: Path,
    fit_views: Path,
    surface_path: Path,
    output: Path,
    outer_fold: int,
    seed: int,
) -> dict[str, Any]:
    """Run the exact production grids for one outer-fold/seed timing unit."""
    if output.exists():
        raise HistologyBaselineError(f"refusing to overwrite fit output: {output}")
    surface = json.loads(surface_path.read_text(encoding="utf-8"))
    if surface.get("production_fit_authorized") is not False:
        raise HistologyBaselineError("validated design no-fit marker changed")
    frozen_seeds = [int(value) for value in surface["resampling"]["model_seed_roster"]]
    if outer_fold not in range(5) or seed not in frozen_seeds:
        raise HistologyBaselineError("preflight fold or seed is outside the frozen roster")
    _, participants = read_tsv(molecular / "participant_axis.tsv")
    _, fold_rows = read_tsv(folds / "participant_outer_folds.tsv")
    participant_ids = [row["participant_id"] for row in participants]
    outer = np.asarray([int(row["outer_fold"]) for row in fold_rows], dtype=int)
    if participant_ids != [row["participant_id"] for row in fold_rows] or len(participant_ids) != 99:
        raise HistologyBaselineError("participant and fold axes differ")
    output.mkdir(parents=True)
    receipts_root = output / "fit_receipts"
    receipts_root.mkdir()
    prediction_root = output / "predictions"
    prediction_root.mkdir()
    result = _fit_outer_bundle(
        outer_fold=outer_fold,
        seeds=[seed],
        molecular=molecular,
        folds=folds,
        fit_views=fit_views,
        surface=surface,
        output=receipts_root,
    )[0]
    indices = np.asarray(result["test_indices"], dtype=int)
    held_ids = [participant_ids[index] for index in indices]
    held_outer = outer[indices]
    for model_id in MODEL_IDS:
        write_tsv(
            prediction_root / f"{model_id}.tsv",
            PREDICTION_FIELDS,
            _prediction_rows(held_ids, held_outer, result["results"][model_id]),
        )
    receipt = {
        "schema_version": "masld-bench-gse267145-histology-baseline-fit-preflight-v1",
        "status": "passed_timing_predictions_unscored",
        "production_predictions_complete": False,
        "exact_production_grids_used": True,
        "outer_fold": outer_fold,
        "seed": seed,
        "training_participants": 99 - len(indices),
        "query_participants": len(indices),
        "logical_outer_seed_tasks": 1,
        "model_family_fits": len(MODEL_IDS),
        "prediction_files": len(MODEL_IDS),
        "fit_counts": result["fit_counts"],
        "elastic_logistic_solver": result["elastic_logistic_solver"],
        "invalid_candidate_count": result["invalid_candidate_count"],
        "invalid_candidate_reasons": result["invalid_candidate_reasons"],
        "one_standard_error_selection": True,
        "outer_test_outcomes_read": False,
        "recorded_sex_read": False,
        "source_stage5_read_or_selected": False,
        "outer_test_metrics_calculated": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def fit_campaign(
    *,
    molecular: Path,
    folds: Path,
    fit_views: Path,
    surface_path: Path,
    output: Path,
    workers: int,
) -> dict[str, Any]:
    if output.exists():
        raise HistologyBaselineError(f"refusing to overwrite fit output: {output}")
    surface = json.loads(surface_path.read_text(encoding="utf-8"))
    if surface.get("production_fit_authorized") is not False:
        raise HistologyBaselineError("validated design no-fit marker changed")
    _, participants = read_tsv(molecular / "participant_axis.tsv")
    _, fold_rows = read_tsv(folds / "participant_outer_folds.tsv")
    participant_ids = [row["participant_id"] for row in participants]
    outer = np.asarray([int(row["outer_fold"]) for row in fold_rows], dtype=int)
    if participant_ids != [row["participant_id"] for row in fold_rows] or len(participant_ids) != 99:
        raise HistologyBaselineError("participant and fold axes differ")
    rna_raw = np.load(molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    h3_raw = np.load(molecular / "h3k27ac_counts.npy", mmap_mode="r", allow_pickle=False)
    if rna_raw.shape != (99, 42163) or h3_raw.shape != (99, 96460):
        raise HistologyBaselineError("molecular fixture shapes differ")
    output.mkdir(parents=True)
    (output / "fit_receipts").mkdir()
    prediction_root = output / "predictions"
    prediction_root.mkdir()
    seeds = [int(value) for value in surface["resampling"]["model_seed_roster"]]
    assembled: dict[tuple[str, int], dict[str, np.ndarray]] = {
        (model, seed): {
            "stage_test": np.full((99, 3), np.nan),
            "nas_test": np.full(99, np.nan),
            "fibrosis_cumulative_test": np.full(99, np.nan),
            "fibrosis_regression_test": np.full(99, np.nan),
            "group_test": np.full((99, 3), np.nan),
        }
        for model in MODEL_IDS
        for seed in seeds
    }
    tasks = [(outer_fold, seed) for outer_fold in range(5) for seed in seeds]
    worker_root = output / "fit_receipts"
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                _fit_outer_bundle,
                outer_fold=outer_fold,
                seeds=seeds,
                molecular=molecular,
                folds=folds,
                fit_views=fit_views,
                surface=surface,
                output=worker_root,
            )
            for outer_fold in range(5)
        ]
        for future in as_completed(futures):
            for result in future.result():
                indices = result["test_indices"]
                for model_id, payload in result["results"].items():
                    for key, values in payload.items():
                        assembled[(model_id, result["seed"])][key][indices] = values
    for (model_id, seed), payload in assembled.items():
        if any(not np.all(np.isfinite(value)) for value in payload.values()):
            raise HistologyBaselineError("outer-test prediction assembly is incomplete")
        write_tsv(
            prediction_root / f"{model_id}--seed-{seed}.tsv",
            PREDICTION_FIELDS,
            _prediction_rows(participant_ids, outer, payload),
        )
    for model_id in MODEL_IDS:
        payload = {
            key: np.mean(
                np.stack([assembled[(model_id, seed)][key] for seed in seeds], axis=0), axis=0
            )
            for key in next(iter(assembled.values()))
        }
        write_tsv(
            prediction_root / f"{model_id}.tsv",
            PREDICTION_FIELDS,
            _prediction_rows(participant_ids, outer, payload),
        )
    receipt = {
        "schema_version": "masld-bench-gse267145-histology-baseline-fit-v1",
        "status": "passed_predictions_unscored",
        "participants": 99,
        "outer_folds": 5,
        "inner_folds": 4,
        "model_ids": list(MODEL_IDS),
        "model_seeds": seeds,
        "logical_outer_seed_tasks": len(tasks),
        "prediction_files": len(MODEL_IDS) * (len(seeds) + 1),
        "one_standard_error_selection": True,
        "outer_test_outcomes_read": False,
        "recorded_sex_read": False,
        "source_stage5_read_or_selected": False,
        "single_cell_methods_used": False,
        "h3_coordinates_inferred": False,
        "published_differential_features_used": False,
        "outer_test_metrics_calculated": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecular", required=True, type=Path)
    parser.add_argument("--molecular-artifacts-sha256", required=True)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--folds-artifacts-sha256", required=True)
    parser.add_argument("--fit-views", required=True, type=Path)
    parser.add_argument("--fit-views-artifacts-sha256", required=True)
    parser.add_argument("--surface", required=True, type=Path)
    parser.add_argument("--surface-sha256", required=True)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--preflight-outer-fold", type=int)
    parser.add_argument("--preflight-seed", type=int)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    for root, expected in (
        (arguments.molecular, arguments.molecular_artifacts_sha256),
        (arguments.folds, arguments.folds_artifacts_sha256),
        (arguments.fit_views, arguments.fit_views_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise HistologyBaselineError(f"input ARTIFACTS SHA differs: {root}")
    if sha256_file(arguments.surface) != arguments.surface_sha256:
        raise HistologyBaselineError("benchmark surface SHA differs")
    preflight_values = (arguments.preflight_outer_fold, arguments.preflight_seed)
    if (preflight_values[0] is None) != (preflight_values[1] is None):
        raise HistologyBaselineError("preflight outer fold and seed must be supplied together")
    if preflight_values[0] is None:
        result = fit_campaign(
            molecular=arguments.molecular,
            folds=arguments.folds,
            fit_views=arguments.fit_views,
            surface_path=arguments.surface,
            output=arguments.output,
            workers=arguments.workers,
        )
    else:
        result = fit_preflight(
            molecular=arguments.molecular,
            folds=arguments.folds,
            fit_views=arguments.fit_views,
            surface_path=arguments.surface,
            output=arguments.output,
            outer_fold=int(preflight_values[0]),
            seed=int(preflight_values[1]),
        )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
