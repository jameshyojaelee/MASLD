#!/usr/bin/env python
"""Regularized donor-held-out RNA-to-ATAC profile baselines.

Preparation delegates to the frozen classical adapter. Fit receives training
RNA and ATAC plus query RNA only through the prepared action. The original
paired HDF5 is withheld from fit and predict by the campaign control plane.
Hyperparameter selection is grouped by donor inside each outer training fold.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence


class RNAATACRegularizedError(RuntimeError):
    """Raised when a regularized RNA-to-ATAC contract is violated."""


def _load_classical() -> Any:
    path = Path(__file__).with_name("rna_atac_classical.py")
    spec = importlib.util.spec_from_file_location(
        "masld_bench_standalone_rna_atac_classical_for_regularized", path
    )
    if spec is None or spec.loader is None:
        raise RNAATACRegularizedError(f"cannot load classical adapter: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = _load_classical()

MODEL_IDS = ("shrunken_pseudobulk", "shuffled_context", "trans_only")
base.MODEL_IDS = MODEL_IDS
LINEAGES = base.LINEAGES
TASK_ID = base.TASK_ID
DATASET_ID = base.DATASET_ID
RUNTIME_ID = base.RUNTIME_ID
RECEIPT_SCHEMA = base.RECEIPT_SCHEMA
PREDICTION_SCHEMA = base.PREDICTION_SCHEMA

INNER_FOLDS = 3
INNER_SPLIT_SEED = 20260823
CONTEXT_PERMUTATION_SEED = 20260824
CONTEXT_RANKS = (4, 8, 16)
RIDGE_LAMBDAS = (0.1, 1.0, 10.0, 100.0)
SHRINKAGE_WEIGHTS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
ATAC_LOG_TARGET_SUM = 10000.0


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row))


def _write_array(path: Path, value: Any) -> None:
    import numpy as np

    with path.open("xb") as handle:
        np.save(handle, value, allow_pickle=False)


def _receipt(
    *,
    action: str,
    request: Mapping[str, Any],
    output: Path,
    artifacts: Sequence[Path],
    environment_sha256: str,
    metadata: Mapping[str, Any],
) -> None:
    base._write_json(
        output / "adapter_receipt.json",
        {
            "schema_version": RECEIPT_SCHEMA,
            "action": action,
            "run_id": request["run_id"],
            "status": "complete",
            "artifacts": [
                base._artifact_record(path, relative_to=output) for path in artifacts
            ],
            "metadata": {
                "adapter": "rna_atac_regularized_v1",
                "runtime_id": RUNTIME_ID,
                "environment_artifact_sha256": environment_sha256,
                "fit_dataset_ids": list(request["fit_dataset_ids"]),
                **dict(metadata),
            },
        },
    )


def _dense_profiles(matrix: Any) -> Any:
    import numpy as np

    values = base._row_normalize(matrix).toarray()
    if (
        values.ndim != 2
        or np.any(~np.isfinite(values))
        or np.any(values < 0)
        or np.any(values.sum(axis=1) <= 0)
    ):
        raise RNAATACRegularizedError("ATAC profile matrix is invalid")
    return values


def _profile_loss(observed: Any, predicted: Any, pseudocount: float) -> float:
    import numpy as np

    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64) + float(pseudocount)
    score /= score.sum()
    if truth.shape != score.shape or np.any(truth < 0) or np.any(score <= 0):
        raise RNAATACRegularizedError("profile loss received invalid values")
    positive = truth > 0
    return float(-np.sum(truth[positive] * np.log(score[positive])))


def _mean_profiles(
    profiles: Any, rows: Sequence[Mapping[str, str]], indices: Sequence[int]
) -> tuple[Any, dict[str, Any]]:
    import numpy as np

    selected = np.asarray(indices, dtype=np.int64)
    if selected.size == 0:
        raise RNAATACRegularizedError("profile mean received an empty training set")
    global_profile = np.asarray(profiles[selected].mean(axis=0), dtype=np.float64)
    lineage_profiles: dict[str, Any] = {}
    for lineage in LINEAGES:
        lineage_indices = [
            index for index in selected if rows[int(index)]["lineage"] == lineage
        ]
        if not lineage_indices:
            raise RNAATACRegularizedError(f"training fold lacks lineage {lineage}")
        lineage_profiles[lineage] = np.asarray(
            profiles[lineage_indices].mean(axis=0), dtype=np.float64
        )
    return global_profile, lineage_profiles


def choose_shrinkage_weight(
    profiles: Any,
    rows: Sequence[Mapping[str, str]],
    *,
    pseudocount: float,
) -> tuple[float, list[dict[str, Any]]]:
    """Select one lineage-to-global weight using grouped inner donor CV."""

    donor_folds = {
        row["donor_id"]: base.fold_index(
            row["donor_id"], seed=INNER_SPLIT_SEED, outer_folds=INNER_FOLDS
        )
        for row in rows
    }
    losses = {weight: 0.0 for weight in SHRINKAGE_WEIGHTS}
    counts = {weight: 0 for weight in SHRINKAGE_WEIGHTS}
    for fold in range(INNER_FOLDS):
        train = [
            index
            for index, row in enumerate(rows)
            if donor_folds[row["donor_id"]] != fold
        ]
        validation = [
            index
            for index, row in enumerate(rows)
            if donor_folds[row["donor_id"]] == fold
        ]
        if not train or not validation:
            raise RNAATACRegularizedError("inner donor split is empty")
        global_profile, lineage_profiles = _mean_profiles(profiles, rows, train)
        for weight in SHRINKAGE_WEIGHTS:
            for index in validation:
                lineage = rows[index]["lineage"]
                prediction = (
                    (1.0 - weight) * global_profile
                    + weight * lineage_profiles[lineage]
                )
                losses[weight] += _profile_loss(
                    profiles[index], prediction, pseudocount
                )
                counts[weight] += 1
    records = [
        {
            "weight": format(weight, ".17g"),
            "mean_cross_entropy": format(losses[weight] / counts[weight], ".17g"),
            "validation_profiles": counts[weight],
        }
        for weight in SHRINKAGE_WEIGHTS
    ]
    selected = min(
        SHRINKAGE_WEIGHTS,
        key=lambda weight: (losses[weight] / counts[weight], weight),
    )
    return selected, records


def fit_shrunken_profiles(
    profiles: Any,
    rows: Sequence[Mapping[str, str]],
    *,
    pseudocount: float,
) -> tuple[Any, float, list[dict[str, Any]]]:
    import numpy as np

    weight, records = choose_shrinkage_weight(
        profiles, rows, pseudocount=pseudocount
    )
    global_profile, lineage_profiles = _mean_profiles(
        profiles, rows, range(len(rows))
    )
    fitted = np.vstack(
        [
            (1.0 - weight) * global_profile + weight * lineage_profiles[lineage]
            for lineage in LINEAGES
        ]
    )
    fitted /= fitted.sum(axis=1, keepdims=True)
    return fitted, weight, records


def context_permutation_indices(
    rows: Sequence[Mapping[str, str]], *, seed: int
) -> tuple[int, ...]:
    """Return a deterministic no-fixed-point donor permutation within lineage."""

    mapping = list(range(len(rows)))
    for lineage in LINEAGES:
        destinations = sorted(
            [index for index, row in enumerate(rows) if row["lineage"] == lineage],
            key=lambda index: (rows[index]["donor_id"], index),
        )
        if len(destinations) < 2:
            raise RNAATACRegularizedError(
                f"cannot permute fewer than two donors in {lineage}"
            )
        offset = 1 + int.from_bytes(
            sha256(f"{seed}\0{lineage}".encode()).digest()[:8], "big"
        ) % (len(destinations) - 1)
        sources = destinations[offset:] + destinations[:offset]
        for destination, source in zip(destinations, sources, strict=True):
            mapping[destination] = source
            if rows[destination]["lineage"] != rows[source]["lineage"]:
                raise RNAATACRegularizedError("context permutation crossed lineages")
    if any(index == source for index, source in enumerate(mapping)):
        raise RNAATACRegularizedError("context permutation retained a fixed point")
    return tuple(mapping)


def _fit_rna_transform(matrix: Any, *, n_hvg: int, max_rank: int) -> dict[str, Any]:
    import numpy as np

    normalized = base._log_normalize(matrix, 10000.0)
    selected = base._select_context_hvgs(
        normalized, min(int(n_hvg), normalized.shape[1])
    )
    values = normalized[:, selected].toarray()
    mean = values.mean(axis=0)
    scale = values.std(axis=0)
    scale[scale == 0] = 1.0
    standardized = (values - mean) / scale
    _, _, right = np.linalg.svd(standardized, full_matrices=False)
    rank = min(int(max_rank), right.shape[0], max(1, standardized.shape[0] - 1))
    components = right[:rank]
    scores = standardized @ components.T
    pc_scale = scores.std(axis=0)
    pc_scale[pc_scale == 0] = 1.0
    scores /= pc_scale
    return {
        "selected": np.asarray(selected, dtype=np.int64),
        "mean": mean,
        "scale": scale,
        "components": components,
        "pc_scale": pc_scale,
        "scores": scores,
    }


def _apply_rna_transform(matrix: Any, state: Mapping[str, Any]) -> Any:
    normalized = base._log_normalize(matrix, 10000.0)
    values = normalized[:, state["selected"]].toarray()
    standardized = (values - state["mean"]) / state["scale"]
    return (standardized @ state["components"].T) / state["pc_scale"]


def _design(scores: Any, rows: Sequence[Mapping[str, str]], rank: int) -> Any:
    import numpy as np

    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}
    one_hot = np.zeros((len(rows), len(LINEAGES)), dtype=np.float64)
    for index, row in enumerate(rows):
        one_hot[index, lineage_index[row["lineage"]]] = 1.0
    return np.hstack((one_hot, np.asarray(scores)[:, :rank]))


def _fit_ridge(design: Any, target: Any, *, rank: int, ridge_lambda: float) -> Any:
    import numpy as np

    gram = design.T @ design
    penalty = np.full(gram.shape[0], 1e-8, dtype=np.float64)
    penalty[len(LINEAGES) : len(LINEAGES) + rank] = float(ridge_lambda)
    gram.flat[:: gram.shape[0] + 1] += penalty
    return np.linalg.solve(gram, design.T @ target)


def _profile_predictions(values: Any, fallback: Any, pseudocount: float) -> Any:
    import numpy as np

    upper = math.log1p(ATAC_LOG_TARGET_SUM)
    scores = np.expm1(np.clip(values, 0.0, upper))
    scores = np.maximum(scores, 0.0)
    totals = scores.sum(axis=1)
    for index in np.flatnonzero(totals <= 0):
        scores[index] = fallback
    scores += float(pseudocount)
    scores /= scores.sum(axis=1, keepdims=True)
    return scores


def _context_grid(
    rna: Any,
    profiles: Any,
    rows: Sequence[Mapping[str, str]],
    *,
    n_hvg: int,
    pseudocount: float,
    shuffled: bool,
) -> tuple[int, float, list[dict[str, Any]]]:
    import numpy as np

    donor_folds = {
        row["donor_id"]: base.fold_index(
            row["donor_id"], seed=INNER_SPLIT_SEED, outer_folds=INNER_FOLDS
        )
        for row in rows
    }
    loss = {(rank, value): 0.0 for rank in CONTEXT_RANKS for value in RIDGE_LAMBDAS}
    count = {(rank, value): 0 for rank in CONTEXT_RANKS for value in RIDGE_LAMBDAS}
    for fold in range(INNER_FOLDS):
        train = np.asarray(
            [
                index
                for index, row in enumerate(rows)
                if donor_folds[row["donor_id"]] != fold
            ],
            dtype=np.int64,
        )
        validation = np.asarray(
            [
                index
                for index, row in enumerate(rows)
                if donor_folds[row["donor_id"]] == fold
            ],
            dtype=np.int64,
        )
        if train.size == 0 or validation.size == 0:
            raise RNAATACRegularizedError("inner context split is empty")
        train_rows = [rows[int(index)] for index in train]
        validation_rows = [rows[int(index)] for index in validation]
        transform = _fit_rna_transform(
            rna[train], n_hvg=n_hvg, max_rank=max(CONTEXT_RANKS)
        )
        train_scores = transform["scores"]
        if shuffled:
            permutation = context_permutation_indices(
                train_rows, seed=CONTEXT_PERMUTATION_SEED + fold
            )
            train_scores = train_scores[list(permutation)]
        validation_scores = _apply_rna_transform(rna[validation], transform)
        target = np.log1p(profiles[train] * ATAC_LOG_TARGET_SUM)
        fallback = profiles[train].mean(axis=0)
        fallback /= fallback.sum()
        for rank in CONTEXT_RANKS:
            train_design = _design(train_scores, train_rows, rank)
            validation_design = _design(validation_scores, validation_rows, rank)
            for ridge_lambda in RIDGE_LAMBDAS:
                coefficients = _fit_ridge(
                    train_design,
                    target,
                    rank=rank,
                    ridge_lambda=ridge_lambda,
                )
                predictions = _profile_predictions(
                    validation_design @ coefficients, fallback, pseudocount
                )
                key = (rank, ridge_lambda)
                for position, observed_index in enumerate(validation):
                    loss[key] += _profile_loss(
                        profiles[int(observed_index)],
                        predictions[position],
                        pseudocount,
                    )
                    count[key] += 1
    records = [
        {
            "pca_rank": rank,
            "ridge_lambda": format(ridge_lambda, ".17g"),
            "mean_cross_entropy": format(
                loss[(rank, ridge_lambda)] / count[(rank, ridge_lambda)], ".17g"
            ),
            "validation_profiles": count[(rank, ridge_lambda)],
        }
        for rank in CONTEXT_RANKS
        for ridge_lambda in RIDGE_LAMBDAS
    ]
    selected_rank, selected_lambda = min(
        loss,
        key=lambda key: (
            loss[key] / count[key],
            key[0],
            -key[1],
        ),
    )
    return selected_rank, selected_lambda, records


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np
    from scipy import sparse

    parameters, model_id = base._validate_request(request, "fit")
    _, environment_sha256 = base._validate_environment(request)
    prepared = base._prior_output(request_path, request, "prepare")
    train_rna = base._load_sparse(prepared / "training_rna.npz")
    train_atac = base._load_sparse(prepared / "training_atac.npz")
    fields, rows = base._read_tsv(prepared / "training_rows.tsv")
    if (
        fields != ("donor_id", "lineage", "n_nuclei")
        or len(rows) != train_rna.shape[0]
        or train_rna.shape[0] != train_atac.shape[0]
    ):
        raise RNAATACRegularizedError("prepared training artifacts differ")
    profiles = _dense_profiles(train_atac)
    pseudocount = float(parameters["profile_pseudocount"])
    artifacts: list[Path] = []
    state: dict[str, Any]
    if model_id == "shrunken_pseudobulk":
        fitted, weight, selection = fit_shrunken_profiles(
            profiles, rows, pseudocount=pseudocount
        )
        profiles_path = output / "profiles.npz"
        base._save_sparse(profiles_path, sparse.csr_matrix(fitted))
        index_path = output / "profile_index.tsv"
        _write_tsv(
            index_path,
            ("profile_id", "lineage", "donor_id"),
            (
                {
                    "profile_id": lineage,
                    "lineage": lineage,
                    "donor_id": "training_shrunken_mean",
                }
                for lineage in LINEAGES
            ),
        )
        selection_path = output / "shrinkage_selection.tsv"
        _write_tsv(
            selection_path,
            ("weight", "mean_cross_entropy", "validation_profiles"),
            selection,
        )
        artifacts.extend((profiles_path, index_path, selection_path))
        state = {
            "model_type": "training_only_cv_shrunken_lineage_profile",
            "selected_shrinkage_weight": weight,
            "profiles": base._artifact_record(profiles_path, relative_to=output),
            "profile_index": base._artifact_record(index_path, relative_to=output),
            "selection": base._artifact_record(selection_path, relative_to=output),
        }
    else:
        shuffled = model_id == "shuffled_context"
        rank, ridge_lambda, selection = _context_grid(
            train_rna,
            profiles,
            rows,
            n_hvg=int(parameters["n_context_hvg"]),
            pseudocount=pseudocount,
            shuffled=shuffled,
        )
        transform = _fit_rna_transform(
            train_rna, n_hvg=int(parameters["n_context_hvg"]), max_rank=rank
        )
        fit_scores = transform["scores"]
        permutation: tuple[int, ...] | None = None
        if shuffled:
            permutation = context_permutation_indices(
                rows, seed=CONTEXT_PERMUTATION_SEED
            )
            fit_scores = fit_scores[list(permutation)]
        design = _design(fit_scores, rows, rank)
        target = np.log1p(profiles * ATAC_LOG_TARGET_SUM)
        coefficients = _fit_ridge(
            design,
            target,
            rank=rank,
            ridge_lambda=ridge_lambda,
        )
        fallback = profiles.mean(axis=0)
        fallback /= fallback.sum()
        arrays = {
            "selected_gene_indices.npy": transform["selected"],
            "gene_mean.npy": transform["mean"],
            "gene_scale.npy": transform["scale"],
            "components.npy": transform["components"],
            "pc_scale.npy": transform["pc_scale"],
            "coefficients.npy": coefficients,
            "fallback_profile.npy": fallback,
        }
        array_records = {}
        for name, value in arrays.items():
            path = output / name
            _write_array(path, value)
            artifacts.append(path)
            array_records[name] = base._artifact_record(path, relative_to=output)
        selection_path = output / "context_selection.tsv"
        _write_tsv(
            selection_path,
            ("pca_rank", "ridge_lambda", "mean_cross_entropy", "validation_profiles"),
            selection,
        )
        artifacts.append(selection_path)
        state = {
            "model_type": "donor_grouped_nested_cv_rna_pca_ridge",
            "selected_pca_rank": rank,
            "selected_ridge_lambda": ridge_lambda,
            "training_context_permuted_within_lineage": shuffled,
            "context_permutation_sha256": (
                base._canonical_hash(list(permutation))
                if permutation is not None
                else None
            ),
            "arrays": array_records,
            "selection": base._artifact_record(selection_path, relative_to=output),
        }
    model_path = output / "fitted_model.json"
    base._write_json(
        model_path,
        {
            "schema_version": "masld-bench-rna-atac-regularized-model-v1",
            "run_id": request["run_id"],
            "model_id": model_id,
            "held_out_fold": int(request["run_spec"]["fold"]),
            "training_profile_count": len(rows),
            "training_donor_count": len({row["donor_id"] for row in rows}),
            "peak_count": train_atac.shape[1],
            "parameter_sha256": base._canonical_hash(dict(parameters)),
            "environment_lock_sha256": environment_sha256,
            "inner_folds": INNER_FOLDS,
            "inner_split_seed": INNER_SPLIT_SEED,
            "query_rna_used_for_fit": False,
            "held_atac_used_for_fit": False,
            **state,
        },
    )
    artifacts.append(model_path)
    _receipt(
        action="fit",
        request=request,
        output=output,
        artifacts=artifacts,
        environment_sha256=environment_sha256,
        metadata={
            "model_id": model_id,
            "held_out_fold": int(request["run_spec"]["fold"]),
            "query_rna_used_for_fit": False,
            "held_atac_used_for_fit": False,
        },
    )


def _load_model(fitted: Path, request: Mapping[str, Any], model_id: str) -> Mapping[str, Any]:
    try:
        model = json.loads((fitted / "fitted_model.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RNAATACRegularizedError("fitted model manifest is invalid") from error
    if (
        not isinstance(model, Mapping)
        or model.get("schema_version") != "masld-bench-rna-atac-regularized-model-v1"
        or model.get("run_id") != request["run_id"]
        or model.get("model_id") != model_id
        or model.get("parameter_sha256")
        != base._canonical_hash(dict(request["run_spec"]["hyperparameters"]))
        or model.get("held_atac_used_for_fit") is not False
    ):
        raise RNAATACRegularizedError("fitted model binding differs")
    return model


def _context_profiles(
    *,
    fitted: Path,
    prepared: Path,
    model: Mapping[str, Any],
    query_rows: Sequence[Mapping[str, str]],
    pseudocount: float,
) -> Any:
    arrays = model.get("arrays")
    if not isinstance(arrays, Mapping):
        raise RNAATACRegularizedError("context model arrays are missing")
    transform = {
        "selected": base._load_array(fitted, arrays["selected_gene_indices.npy"]).astype("int64"),
        "mean": base._load_array(fitted, arrays["gene_mean.npy"]),
        "scale": base._load_array(fitted, arrays["gene_scale.npy"]),
        "components": base._load_array(fitted, arrays["components.npy"]),
        "pc_scale": base._load_array(fitted, arrays["pc_scale.npy"]),
    }
    query_rna = base._load_sparse(prepared / "query_rna.npz")
    scores = _apply_rna_transform(query_rna, transform)
    rank = int(model["selected_pca_rank"])
    design = _design(scores, query_rows, rank)
    coefficients = base._load_array(fitted, arrays["coefficients.npy"])
    fallback = base._load_array(fitted, arrays["fallback_profile.npy"])
    return _profile_predictions(design @ coefficients, fallback, pseudocount)


def predict(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np

    parameters, model_id = base._validate_request(request, "predict")
    _, environment_sha256 = base._validate_environment(request)
    prepared = base._prior_output(request_path, request, "prepare")
    fitted = base._prior_output(request_path, request, "fit")
    model = _load_model(fitted, request, model_id)
    query_fields, query_rows = base._read_tsv(prepared / "query_rows.tsv")
    peak_fields, peak_rows = base._read_tsv(prepared / "selected_peaks.tsv")
    if query_fields != ("donor_id", "lineage", "n_nuclei"):
        raise RNAATACRegularizedError("query row schema differs")
    if peak_fields != (
        "selected_index",
        "source_index",
        "peak_id",
        "chromosome",
        "bed_start",
        "bed_end",
    ):
        raise RNAATACRegularizedError("selected peak schema differs")
    pseudocount = float(parameters["profile_pseudocount"])
    if model_id == "shrunken_pseudobulk":
        profiles = base._load_sparse(fitted / "profiles.npz").toarray()
        fields, profile_rows = base._read_tsv(fitted / "profile_index.tsv")
        if fields != ("profile_id", "lineage", "donor_id"):
            raise RNAATACRegularizedError("shrunken profile index differs")
        by_lineage = {row["lineage"]: index for index, row in enumerate(profile_rows)}
        if set(by_lineage) != set(LINEAGES):
            raise RNAATACRegularizedError("shrunken lineage roster differs")
        predicted_profiles = np.vstack(
            [profiles[by_lineage[row["lineage"]]] for row in query_rows]
        )
    else:
        predicted_profiles = _context_profiles(
            fitted=fitted,
            prepared=prepared,
            model=model,
            query_rows=query_rows,
            pseudocount=pseudocount,
        )
    if predicted_profiles.shape != (len(query_rows), len(peak_rows)):
        raise RNAATACRegularizedError("predicted profile axis differs")
    namespace = str(parameters["join_namespace"])
    prediction_fields = (
        "row_hash",
        "donor_hash",
        "block_hash",
        "stratum",
        "predicted",
    )
    prediction_rows: list[dict[str, str]] = []
    row_id_rows: list[dict[str, str]] = []
    for query_index, query in enumerate(query_rows):
        donor_hash = base.join_hash(namespace, "unit", query["donor_id"])
        for peak_index, peak in enumerate(peak_rows):
            row_hash = base.join_hash(
                namespace,
                "row",
                f"{query['donor_id']}\0{query['lineage']}\0{peak['peak_id']}",
            )
            prediction_rows.append(
                {
                    "row_hash": row_hash,
                    "donor_hash": donor_hash,
                    "block_hash": base.join_hash(
                        namespace, "block", peak["chromosome"]
                    ),
                    "stratum": query["lineage"],
                    "predicted": format(
                        float(predicted_profiles[query_index, peak_index])
                        + pseudocount,
                        ".17g",
                    ),
                }
            )
            row_id_rows.append({"row_hash": row_hash, "donor_hash": donor_hash})
    prediction_rows.sort(key=lambda row: row["row_hash"])
    row_id_rows.sort(key=lambda row: row["row_hash"])
    prediction_path = output / "predictions.tsv"
    row_ids_path = output / "row_ids.tsv"
    _write_tsv(prediction_path, prediction_fields, prediction_rows)
    _write_tsv(row_ids_path, ("row_hash", "donor_hash"), row_id_rows)
    standardized = {
        **base._artifact_record(prediction_path, relative_to=output),
        "media_type": "text/tab-separated-values",
        "role": f"standardized_prediction_table:{TASK_ID}",
    }
    row_artifact = {
        **base._artifact_record(row_ids_path, relative_to=output),
        "media_type": "text/tab-separated-values",
        "role": f"prediction_row_ids:{TASK_ID}",
    }
    source_join = base._canonical_hash(
        {
            "task_id": TASK_ID,
            "dataset_ids": [DATASET_ID],
            "split_id": "donor_outer",
            "row_id_field": "row_hash",
            "unit_id_field": "donor_hash",
            "unit_id_namespace": namespace,
            "biological_unit": "donor",
        }
    )
    bundle = {
        "schema_version": PREDICTION_SCHEMA,
        "bundle_id": f"{model_id}-{str(request['run_id'])[:16]}",
        "run_id": request["run_id"],
        "task_id": TASK_ID,
        "model_id": model_id,
        "dataset_ids": [DATASET_ID],
        "split_id": "donor_outer",
        "artifacts": [standardized, row_artifact],
        "standardized_table": standardized,
        "row_ids": row_artifact,
        "n_predictions": len(prediction_rows),
        "row_id_field": "row_hash",
        "unit_id_field": "donor_hash",
        "unit_id_namespace": namespace,
        "biological_unit": "donor",
        "table_schema_sha256": base._canonical_hash(
            {"format": "tsv", "fields": list(prediction_fields)}
        ),
        "source_join_key_sha256": source_join,
        "format_version": "tsv-v1",
        "missing_state": "observed",
        "metadata": {
            "held_out_fold": int(request["run_spec"]["fold"]),
            "query_stratum_count": len(query_rows),
            "selected_peak_count": len(peak_rows),
            "prediction_scale": "positive_depth_free_rna_conditioned_profile",
            "profile_pseudocount": pseudocount,
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
            "smoke_only": True,
        },
    }
    bundle_path = output / "prediction_bundle.json"
    base._write_json(bundle_path, bundle)
    _receipt(
        action="predict",
        request=request,
        output=output,
        artifacts=(prediction_path, row_ids_path, bundle_path),
        environment_sha256=environment_sha256,
        metadata={
            "model_id": model_id,
            "held_out_fold": int(request["run_spec"]["fold"]),
            "prediction_row_count": len(prediction_rows),
            "prediction_donor_count": len(
                {row["donor_hash"] for row in prediction_rows}
            ),
            "held_atac_input_exposed": False,
            "observed_atac_exported": False,
        },
    )


def prepare(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    base.prepare(request_path, request, output)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", required=True, choices=("prepare", "fit", "predict"))
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    try:
        request = json.loads(arguments.request.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RNAATACRegularizedError("adapter request is invalid JSON") from error
    if not isinstance(request, Mapping):
        raise RNAATACRegularizedError("adapter request must be an object")
    if arguments.output.exists():
        raise RNAATACRegularizedError(
            f"adapter output already exists: {arguments.output}"
        )
    arguments.output.mkdir(parents=True, exist_ok=False)
    if arguments.action == "prepare":
        prepare(arguments.request, request, arguments.output)
    elif arguments.action == "fit":
        fit(arguments.request, request, arguments.output)
    else:
        predict(arguments.request, request, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
