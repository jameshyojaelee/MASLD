#!/usr/bin/env python3
"""Diagnose elastic-net optimization on one frozen study-held Atlas fold."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import time
import warnings

import numpy as np

from scripts.fit_predict_cell_baselines_study_50000 import (
    N_PCA,
    N_TOP_HVG,
    ROSTER,
    TARGET_SUM,
    _log_normalize,
    _select_hvgs,
    donor_class_weights,
    standardize_projection,
)


CONFIGS = (
    {"config_id": "float64_tol1e-4_c1", "dtype": "float64", "tol": 1.0e-4, "c": 1.0, "max_iter": 10_000},
    {"config_id": "float64_tol1e-3_c1", "dtype": "float64", "tol": 1.0e-3, "c": 1.0, "max_iter": 10_000},
    {"config_id": "float64_tol1e-2_c1", "dtype": "float64", "tol": 1.0e-2, "c": 1.0, "max_iter": 10_000},
    {"config_id": "float64_tol1e-3_c0.1", "dtype": "float64", "tol": 1.0e-3, "c": 0.1, "max_iter": 10_000},
    {"config_id": "float32_tol1e-3_c1", "dtype": "float32", "tol": 1.0e-3, "c": 1.0, "max_iter": 10_000},
    {"config_id": "float64_tol1e-3_c1_max50000", "dtype": "float64", "tol": 1.0e-3, "c": 1.0, "max_iter": 50_000},
)


class ElasticDiagnosticError(RuntimeError):
    """Raised when the frozen diagnostic input or fit differs."""


def select_configs(config_ids: list[str] | None) -> tuple[dict[str, object], ...]:
    lookup = {str(config["config_id"]): config for config in CONFIGS}
    if not config_ids:
        return tuple(CONFIGS[:5])
    if len(config_ids) != len(set(config_ids)) or any(value not in lookup for value in config_ids):
        raise ElasticDiagnosticError("diagnostic configuration roster differs")
    return tuple(lookup[value] for value in config_ids)


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def quantiles(values: np.ndarray) -> dict[str, float]:
    observed = np.asarray(values, dtype=np.float64)
    return {
        label: float(np.quantile(observed, probability))
        for label, probability in (
            ("min", 0.0),
            ("q01", 0.01),
            ("q50", 0.5),
            ("q99", 0.99),
            ("max", 1.0),
        )
    }


def main() -> int:
    import anndata
    from sklearn.decomposition import PCA
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--split", required=True, type=Path)
    parser.add_argument("--split-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--outer-fold", type=int, default=0, choices=range(5))
    parser.add_argument("--config-id", action="append", dest="config_ids")
    arguments = parser.parse_args()
    if (
        sha256_file(arguments.source / "ARTIFACTS.json")
        != arguments.source_artifacts_sha256
        or sha256_file(arguments.split / "ARTIFACTS.json")
        != arguments.split_artifacts_sha256
    ):
        raise ElasticDiagnosticError("source or split ARTIFACTS SHA-256 differs")
    split_rows = read_tsv(arguments.split / "row_outer_folds.tsv")
    adata = anndata.read_h5ad(
        arguments.source / "resource_atlas_frozen_screen_50000.h5ad"
    )
    if len(split_rows) != 50_000 or adata.n_obs != 50_000:
        raise ElasticDiagnosticError("frozen row count differs")
    row_ids = np.asarray(list(map(str, adata.obs["row_id"])), dtype=str)
    if row_ids.tolist() != [row["row_id"] for row in split_rows]:
        raise ElasticDiagnosticError("split and matrix row order differ")
    outer = np.asarray([int(row["outer_fold"]) for row in split_rows])
    donors = np.asarray([row["donor_id"] for row in split_rows], dtype=str)
    labels = np.asarray(list(map(str, adata.obs["broad_label"])), dtype=str)
    train = np.flatnonzero(outer != arguments.outer_fold)
    if set(labels[train]) != set(ROSTER):
        raise ElasticDiagnosticError("outer-fold training class roster differs")
    normalized = _log_normalize(adata.X[train], TARGET_SUM)
    gene_ids = list(map(str, adata.var["ensembl_id"]))
    selected, _ = _select_hvgs(normalized, gene_ids, n_top=N_TOP_HVG, n_bins=20)
    dense = normalized[:, selected].toarray().astype(np.float32, copy=False)
    pca = PCA(
        n_components=N_PCA,
        svd_solver="randomized",
        whiten=False,
        random_state=20260824,
    )
    projected = pca.fit_transform(dense).astype(np.float32, copy=False)
    standardized, _, mean, scale = standardize_projection(projected, projected[:1])
    label_to_index = {label: index for index, label in enumerate(ROSTER)}
    target = np.asarray([label_to_index[label] for label in labels[train]])
    weights = np.asarray(
        donor_class_weights(donors[train].tolist(), labels[train].tolist()),
        dtype=np.float64,
    )
    fits: list[dict[str, object]] = []
    train_probabilities: dict[str, np.ndarray] = {}
    for config in select_configs(arguments.config_ids):
        values = standardized.astype(str(config["dtype"]), copy=True)
        classifier = LogisticRegression(
            C=float(config["c"]),
            penalty="elasticnet",
            solver="saga",
            l1_ratio=0.5,
            tol=float(config["tol"]),
            max_iter=int(config["max_iter"]),
            fit_intercept=True,
            random_state=20260824,
        )
        started = time.monotonic()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            classifier.fit(values, target, sample_weight=weights)
        elapsed = time.monotonic() - started
        messages = [
            str(item.message)
            for item in caught
            if issubclass(item.category, ConvergenceWarning)
        ]
        iterations = [int(value) for value in classifier.n_iter_]
        probabilities = np.asarray(classifier.predict_proba(values), dtype=np.float64)
        train_probabilities[str(config["config_id"])] = probabilities
        fits.append(
            {
                **config,
                "iterations": iterations,
                "converged_before_ceiling": max(iterations) < classifier.max_iter,
                "convergence_warnings": messages,
                "elapsed_seconds": elapsed,
                "classes": [int(value) for value in classifier.classes_],
                "coefficients_finite": bool(np.all(np.isfinite(classifier.coef_))),
                "coefficient_zero_fraction": float(np.mean(classifier.coef_ == 0.0)),
                "coefficient_absolute_max": float(np.max(np.abs(classifier.coef_))),
                "training_probability_sha256": sha256(
                    np.ascontiguousarray(probabilities).tobytes()
                ).hexdigest(),
            }
        )
    pairwise_probability_differences: list[dict[str, object]] = []
    roster = sorted(train_probabilities)
    for left_index, left_id in enumerate(roster):
        for right_id in roster[left_index + 1 :]:
            difference = np.abs(train_probabilities[left_id] - train_probabilities[right_id])
            pairwise_probability_differences.append(
                {
                    "left_config_id": left_id,
                    "right_config_id": right_id,
                    "mean_absolute_difference": float(np.mean(difference)),
                    "max_absolute_difference": float(np.max(difference)),
                }
            )
    result = {
        "schema_version": "masld-bench-elastic-net-convergence-diagnostic-v1",
        "status": "completed_diagnostic_no_model_selection",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "outer_fold": arguments.outer_fold,
        "seed": 20260824,
        "training_rows": len(train),
        "training_donors": len(set(donors[train])),
        "training_weight_quantiles": quantiles(weights),
        "projection_scale_quantiles": quantiles(scale),
        "projection_mean_absolute_max": float(np.max(np.abs(mean))),
        "fits": fits,
        "pairwise_training_probability_differences": pairwise_probability_differences,
        "metrics_calculated": False,
        "test_partition_predicted": False,
        "sealed_outcomes_read": False,
    }
    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "diagnostic.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
