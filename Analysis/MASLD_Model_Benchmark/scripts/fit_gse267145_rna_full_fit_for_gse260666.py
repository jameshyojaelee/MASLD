#!/usr/bin/env python3
"""Fit the frozen RNA stage3 source models on all 99 GSE267145 participants."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence
import warnings

import numpy as np


CLASSES = ("NOR", "NAFL", "NASH")


class SourceFullFitError(RuntimeError):
    """Raised when the source full fit would violate its frozen recipe."""


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
            raise SourceFullFitError(f"TSV has no header: {path}")
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


def log1p_cpm(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2 or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise SourceFullFitError("source RNA matrix is not finite and nonnegative")
    totals = values.sum(axis=1)
    if np.any(totals <= 0):
        raise SourceFullFitError("source participant RNA library is empty")
    return np.log1p(values * (1_000_000.0 / totals[:, None]))


def class_weights(labels: np.ndarray) -> dict[int, float]:
    counts = Counter(int(value) for value in labels)
    if set(counts) != {0, 1, 2}:
        raise SourceFullFitError("source fitting partition lacks a stage3 class")
    return {key: len(labels) / (3.0 * count) for key, count in counts.items()}


def select_feature_indices(
    fitting_matrix: np.ndarray,
    stable_gene_ids: Sequence[str],
    feature_count: int,
) -> np.ndarray:
    values = np.asarray(fitting_matrix, dtype=np.float64)
    ids = np.asarray(stable_gene_ids, dtype=str)
    if values.ndim != 2 or values.shape[1] != len(ids):
        raise SourceFullFitError("source feature-selection axis differs")
    minimum_positive = max(3, math.ceil(0.10 * values.shape[0]))
    eligible = np.flatnonzero(np.sum(values > 0, axis=0) >= minimum_positive)
    if len(eligible) < feature_count:
        raise SourceFullFitError("source common axis has too few eligible genes")
    variances = np.var(values[:, eligible], axis=0, ddof=0)
    order = np.lexsort((ids[eligible], -variances))
    return np.asarray(eligible[order[:feature_count]], dtype=np.int64)


def fit_representation(
    matrix: np.ndarray,
    fitting_indices: np.ndarray,
    evaluation_indices: np.ndarray,
    stable_gene_ids: Sequence[str],
    *,
    feature_count: int,
    pca_components: int,
    random_state: int,
) -> dict[str, Any]:
    from sklearn.decomposition import PCA

    fitting = np.asarray(fitting_indices, dtype=np.int64)
    evaluation = np.asarray(evaluation_indices, dtype=np.int64)
    selected = select_feature_indices(matrix[fitting], stable_gene_ids, feature_count)
    fit_dense = np.asarray(matrix[fitting][:, selected], dtype=np.float64)
    center = fit_dense.mean(axis=0)
    scale = fit_dense.std(axis=0)
    scale[scale <= np.finfo(np.float64).eps] = 1.0
    fit_scaled = (fit_dense - center) / scale
    if pca_components >= min(fit_scaled.shape):
        raise SourceFullFitError("source partition cannot support locked PCA size")
    pca = PCA(
        n_components=pca_components,
        svd_solver="randomized",
        whiten=False,
        random_state=random_state,
    )
    fit_projection = pca.fit_transform(fit_scaled)
    if evaluation.size:
        evaluation_dense = np.asarray(matrix[evaluation][:, selected], dtype=np.float64)
        evaluation_projection = pca.transform((evaluation_dense - center) / scale)
    else:
        evaluation_projection = np.empty((0, pca_components), dtype=np.float64)
    return {
        "selected_indices": selected,
        "center": center,
        "scale": scale,
        "pca_mean": np.asarray(pca.mean_, dtype=np.float64),
        "pca_components": np.asarray(pca.components_, dtype=np.float64),
        "fit_projection": np.asarray(fit_projection, dtype=np.float64),
        "evaluation_projection": evaluation_projection,
    }


def fit_linear_svm(x: np.ndarray, y: np.ndarray, *, c_value: float, seed: int):
    from sklearn.svm import LinearSVC

    model = LinearSVC(
        C=c_value,
        dual="auto",
        class_weight=class_weights(y),
        max_iter=20_000,
        random_state=seed,
    )
    model.fit(x, y)
    if int(model.n_iter_) >= 20_000 or tuple(int(value) for value in model.classes_) != (0, 1, 2):
        raise SourceFullFitError("source full-fit linear SVM did not converge")
    return model


def fit_calibrator(raw: np.ndarray, labels: np.ndarray, *, seed: int):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model = LogisticRegression(
            C=1.0,
            l1_ratio=0.0,
            max_iter=5_000,
            class_weight=class_weights(labels),
            random_state=seed,
        )
        model.fit(raw, labels)
    if (
        any(isinstance(item.message, ConvergenceWarning) for item in caught)
        or tuple(int(value) for value in model.classes_) != (0, 1, 2)
    ):
        raise SourceFullFitError("source full-fit SVM calibrator did not converge")
    return model


def fit_elastic_net(
    x: np.ndarray,
    labels: np.ndarray,
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
            max_iter=500_000,
            tol=1e-4,
            class_weight=class_weights(labels),
            random_state=seed,
        )
        model.fit(x, labels)
    if (
        any(isinstance(item.message, ConvergenceWarning) for item in caught)
        or int(np.max(model.n_iter_)) >= 500_000
        or tuple(int(value) for value in model.classes_) != (0, 1, 2)
    ):
        raise SourceFullFitError("source full-fit elastic net did not converge")
    return model


def fit_source_models(
    *,
    source_matrix: np.ndarray,
    stable_gene_ids: Sequence[str],
    labels: np.ndarray,
    outer_folds: np.ndarray,
    seeds: Sequence[int],
    feature_count: int,
    pca_components: int,
    svm_c: float,
    elastic_c: float,
    elastic_l1_ratio: float,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise SourceFullFitError(f"refusing to overwrite source full fit: {output}")
    matrix = np.asarray(source_matrix, dtype=np.float64)
    if matrix.shape[0] != len(labels) or matrix.shape[0] != len(outer_folds):
        raise SourceFullFitError("source participant axes differ")
    if set(int(value) for value in outer_folds) != set(range(5)):
        raise SourceFullFitError("source outer-fold roster differs")
    output.mkdir(parents=True)
    seed_root = output / "seed_models"
    seed_root.mkdir()
    full_indices = np.arange(matrix.shape[0], dtype=np.int64)
    final_feature_ids: list[str] | None = None
    seed_rows: list[dict[str, Any]] = []
    for seed in seeds:
        raw_oof = np.full((matrix.shape[0], 3), np.nan, dtype=np.float64)
        fold_feature_hashes: dict[str, str] = {}
        for outer_fold in range(5):
            fitting = np.flatnonzero(outer_folds != outer_fold)
            evaluation = np.flatnonzero(outer_folds == outer_fold)
            representation = fit_representation(
                matrix,
                fitting,
                evaluation,
                stable_gene_ids,
                feature_count=feature_count,
                pca_components=pca_components,
                random_state=int(seed) + 10_000 * outer_fold,
            )
            svm = fit_linear_svm(
                representation["fit_projection"],
                labels[fitting],
                c_value=svm_c,
                seed=int(seed) + 10_000 * outer_fold + 1_000,
            )
            raw_oof[evaluation] = svm.decision_function(
                representation["evaluation_projection"]
            )
            fold_ids = [
                str(stable_gene_ids[index])
                for index in representation["selected_indices"]
            ]
            fold_feature_hashes[str(outer_fold)] = hashlib.sha256(
                "\n".join(fold_ids).encode("utf-8")
            ).hexdigest()
        if not np.all(np.isfinite(raw_oof)):
            raise SourceFullFitError("source OOF SVM decisions are incomplete")
        calibrator = fit_calibrator(raw_oof, labels, seed=int(seed) + 92_000)
        final_representation = fit_representation(
            matrix,
            full_indices,
            np.asarray([], dtype=np.int64),
            stable_gene_ids,
            feature_count=feature_count,
            pca_components=pca_components,
            random_state=int(seed) + 90_000,
        )
        final_ids = [
            str(stable_gene_ids[index])
            for index in final_representation["selected_indices"]
        ]
        if final_feature_ids is None:
            final_feature_ids = final_ids
        elif final_feature_ids != final_ids:
            raise SourceFullFitError("seed-specific source feature axes differ")
        projection = final_representation["fit_projection"]
        svm = fit_linear_svm(
            projection,
            labels,
            c_value=svm_c,
            seed=int(seed) + 91_000,
        )
        elastic = fit_elastic_net(
            projection,
            labels,
            c_value=elastic_c,
            l1_ratio=elastic_l1_ratio,
            seed=int(seed) + 91_100,
        )
        centroids = np.vstack(
            [projection[labels == class_id].mean(axis=0) for class_id in range(3)]
        )
        current = seed_root / f"seed_{seed}"
        current.mkdir()
        np.savez_compressed(
            current / "model_state.npz",
            center=final_representation["center"],
            scale=final_representation["scale"],
            pca_mean=final_representation["pca_mean"],
            pca_components=final_representation["pca_components"],
            svm_coef=np.asarray(svm.coef_, dtype=np.float64),
            svm_intercept=np.asarray(svm.intercept_, dtype=np.float64),
            svm_calibrator_coef=np.asarray(calibrator.coef_, dtype=np.float64),
            svm_calibrator_intercept=np.asarray(calibrator.intercept_, dtype=np.float64),
            elastic_coef=np.asarray(elastic.coef_, dtype=np.float64),
            elastic_intercept=np.asarray(elastic.intercept_, dtype=np.float64),
            centroid_centroids=np.asarray(centroids, dtype=np.float64),
        )
        state_sha256 = sha256_file(current / "model_state.npz")
        seed_receipt = {
            "schema_version": "masld-bench-gse267145-rna-full-fit-seed-v1",
            "model_seed": int(seed),
            "participants": int(matrix.shape[0]),
            "selected_features": len(final_ids),
            "pca_components": pca_components,
            "svm_c": svm_c,
            "elastic_c": elastic_c,
            "elastic_l1_ratio": elastic_l1_ratio,
            "outer_fold_oof_calibration": True,
            "outer_fold_feature_axis_sha256": fold_feature_hashes,
            "full_fit_feature_axis_sha256": hashlib.sha256(
                "\n".join(final_ids).encode("utf-8")
            ).hexdigest(),
            "model_state_sha256": state_sha256,
            "external_expression_values_read": False,
            "external_labels_read": False,
        }
        with (current / "receipt.json").open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(seed_receipt, indent=2, sort_keys=True) + "\n")
        seed_rows.append(
            {
                "model_seed": seed,
                "model_state_path": f"seed_models/seed_{seed}/model_state.npz",
                "model_state_sha256": state_sha256,
                "selected_feature_axis_sha256": seed_receipt[
                    "full_fit_feature_axis_sha256"
                ],
            }
        )
    assert final_feature_ids is not None
    write_tsv(
        output / "selected_feature_axis.tsv",
        ("feature_index", "stable_gene_id"),
        [
            {"feature_index": index, "stable_gene_id": stable_id}
            for index, stable_id in enumerate(final_feature_ids)
        ],
    )
    write_tsv(
        output / "seed_model_index.tsv",
        tuple(seed_rows[0]),
        seed_rows,
    )
    priors = np.bincount(labels, minlength=3).astype(np.float64)
    priors /= priors.sum()
    receipt = {
        "schema_version": "masld-bench-gse267145-rna-full-fit-v1",
        "status": "passed_source_full99_fit_preprocessing_locked",
        "participants": int(matrix.shape[0]),
        "stage3_class_counts": {
            CLASSES[index]: int(np.sum(labels == index)) for index in range(3)
        },
        "training_stage_distribution": {
            CLASSES[index]: float(priors[index]) for index in range(3)
        },
        "model_seeds": [int(value) for value in seeds],
        "seeds_are_biological_replicates": False,
        "selected_features": len(final_feature_ids),
        "pca_components": pca_components,
        "source_only_outer_fold_oof_svm_calibration": True,
        "external_expression_values_read": False,
        "external_labels_read": False,
        "secondary_source_endpoints_used": False,
        "source_stage5_used": False,
        "recorded_sex_used": False,
        "external_fit_or_calibration_performed": False,
        "source_full_fit_is_a_new_development_score": False,
        "model_selected_or_repaired_from_external_outcomes": False,
        "project_sealed": False,
        "champion_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "clinical_claim_allowed": False,
    }
    with (output / "fit_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def _check_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path) != expected:
        raise SourceFullFitError(f"{label} SHA-256 differs")


def run_full_fit(
    *, benchmark_root: Path, contract_path: Path, contract_sha256: str, output: Path
) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    _check_hash(contract_path, contract_sha256, "source full-fit contract")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "registered_source_full_fit_pending"
        or contract.get("selected_model_id") != "rna_hvg_pca_linear_svm"
        or contract.get("firewall", {}).get("external_labels_read") is not False
    ):
        raise SourceFullFitError("source full-fit contract differs")
    activation_spec = contract["activation"]
    activation_root = benchmark_root / activation_spec["path"]
    _check_hash(
        activation_root / "ARTIFACTS.json",
        activation_spec["artifacts_sha256"],
        "external activation",
    )
    verify_frozen_tree(activation_root)
    activation = json.loads(
        (activation_root / "activation_receipt.json").read_text(encoding="utf-8")
    )
    if (
        activation.get("status") != activation_spec["required_status"]
        or activation.get("selected_source_model_id") != contract["selected_model_id"]
        or activation.get("external_labels_read") is not False
    ):
        raise SourceFullFitError("external activation gate differs")
    selection_spec = contract["selection"]
    selection_root = benchmark_root / selection_spec["path"]
    _check_hash(
        selection_root / "ARTIFACTS.json",
        selection_spec["artifacts_sha256"],
        "RNA selection",
    )
    verify_frozen_tree(selection_root)
    selection = json.loads(
        (selection_root / "selection_receipt.json").read_text(encoding="utf-8")
    )
    if (
        selection.get("status") != selection_spec["required_status"]
        or selection.get("selected_model_id") != contract["selected_model_id"]
    ):
        raise SourceFullFitError("RNA finalist selection differs")
    inputs = contract["source_inputs"]
    molecular = benchmark_root / inputs["molecular_path"]
    outcomes = benchmark_root / inputs["outcomes_path"]
    folds = benchmark_root / inputs["folds_path"]
    for root, expected, label in (
        (molecular, inputs["molecular_artifacts_sha256"], "source molecular"),
        (outcomes, inputs["outcomes_artifacts_sha256"], "source outcomes"),
        (folds, inputs["folds_artifacts_sha256"], "source folds"),
    ):
        _check_hash(root / "ARTIFACTS.json", expected, label)
        verify_frozen_tree(root)
    _, participants = read_tsv(molecular / "participant_axis.tsv")
    _, source_axis = read_tsv(molecular / "rna_feature_axis.tsv")
    _, endpoint_rows = read_tsv(outcomes / "participant_endpoints.tsv")
    _, fold_rows = read_tsv(folds / "participant_outer_folds.tsv")
    participant_ids = [row["participant_id"] for row in participants]
    if (
        len(participant_ids) != inputs["participants"]
        or participant_ids != [row["participant_id"] for row in endpoint_rows]
        or participant_ids != [row["participant_id"] for row in fold_rows]
    ):
        raise SourceFullFitError("source participant axes differ")
    labels = np.asarray(
        [CLASSES.index(row["stage3"]) for row in endpoint_rows], dtype=np.int64
    )
    class_counts = {
        CLASSES[index]: int(np.sum(labels == index)) for index in range(3)
    }
    if class_counts != inputs["stage3_class_counts"]:
        raise SourceFullFitError("source stage3 class census differs")
    outer = np.asarray([int(row["outer_fold"]) for row in fold_rows], dtype=np.int64)
    source_ids = [row["stable_gene_id"] for row in source_axis]
    common_fields, common_rows = read_tsv(activation_root / "common_stable_gene_axis.tsv")
    if "stable_gene_id" not in common_fields:
        raise SourceFullFitError("activation common-gene axis differs")
    common_ids = [row["stable_gene_id"] for row in common_rows]
    if len(common_ids) != activation_spec["common_stable_genes"]:
        raise SourceFullFitError("activation common-gene census differs")
    source_lookup = {value: index for index, value in enumerate(source_ids)}
    if len(source_lookup) != len(source_ids) or any(value not in source_lookup for value in common_ids):
        raise SourceFullFitError("activation gene is absent from source RNA axis")
    raw = np.load(molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    if raw.shape != (inputs["participants"], len(source_ids)):
        raise SourceFullFitError("source RNA matrix shape differs")
    transformed = log1p_cpm(raw)
    common_matrix = transformed[:, [source_lookup[value] for value in common_ids]]
    recipe = contract["recipe"]
    return fit_source_models(
        source_matrix=common_matrix,
        stable_gene_ids=common_ids,
        labels=labels,
        outer_folds=outer,
        seeds=[int(value) for value in recipe["model_seeds"]],
        feature_count=int(recipe["feature_count"]),
        pca_components=int(recipe["pca_components"]),
        svm_c=float(recipe["rna_hvg_pca_linear_svm"]["c"]),
        elastic_c=float(recipe["rna_hvg_pca_elastic_net"]["c"]),
        elastic_l1_ratio=float(
            recipe["rna_hvg_pca_elastic_net"]["l1_ratio"]
        ),
        output=output,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run_full_fit(
        benchmark_root=arguments.benchmark_root,
        contract_path=arguments.contract,
        contract_sha256=arguments.contract_sha256,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
