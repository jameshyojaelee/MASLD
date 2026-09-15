#!/usr/bin/env python3
"""Fit prespecified five-seed GSE281364 heads without calculating metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.linear_model import Ridge

from masld_bench.artifacts import verify_frozen_tree
from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.freeze_gse281364_open_sequence_taskspec import (
    BASELINES,
    CONTEXTS,
    HEADS,
    MISSING_BASELINES,
    MODELS,
    PREDICTION_FIELDS,
    SEEDS,
    TaskSpecError,
    file_sha256,
    load_config,
    safe_path,
    validate_config,
)
from scripts.gse281364_dna_lm_native_contract import apply_projection, head_features


SCHEMA = "masld-bench-gse281364-open-sequence-head-fit-v1"
ROW_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
)
PROJECTION_WIDTH = 256
STD_FLOOR = 1.0e-6
EIGENVALUE_ABSOLUTE_FLOOR = 1.0e-8
EIGENVALUE_RELATIVE_FLOOR = 1.0e-6


class HeadFitError(RuntimeError):
    """Raised when a fit could see held blocks or leave the frozen requirements."""


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise HeadFitError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise HeadFitError(f"{label} must be a JSON object")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        raise HeadFitError(f"cannot read TSV {path}: {error}") from error


def array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(str(array.shape).encode("ascii"))
    digest.update(array.tobytes())
    return digest.hexdigest()


def verify_taskspec(path: Path, expected_sha256: str, config_sha256: str) -> dict[str, Any]:
    if file_sha256(path / "ARTIFACTS.json") != expected_sha256:
        raise HeadFitError("TaskSpec ARTIFACTS hash differs")
    verify_frozen_tree(path)
    receipt = load_json(path / "receipt.json", label="TaskSpec receipt")
    spec = load_json(path / "task_spec.json", label="TaskSpec")
    if (
        receipt.get("status") != "pass_prespecified_taskspec"
        or receipt.get("candidate_count") != 12
        or receipt.get("fixed_seeds") != list(SEEDS)
        or receipt.get("config_sha256") != config_sha256
        or receipt.get("outcome_values_read")
        or receipt.get("model_fit")
        or spec.get("status") != "prespecified_exposed_development_head_campaign"
    ):
        raise HeadFitError("TaskSpec receipt differs")
    return spec


def load_row_universe(
    path: Path,
) -> tuple[list[dict[str, str]], list[str], dict[str, dict[str, str]]]:
    rows = read_tsv(path)
    if not rows or tuple(rows[0]) != ROW_FIELDS or len(rows) != 10330:
        raise HeadFitError("row-universe schema or denominator differs")
    metadata: dict[str, dict[str, str]] = {}
    identities: set[tuple[int, str]] = set()
    seed_counts = {seed: 0 for seed in SEEDS}
    for row in rows:
        seed = int(row["seed"])
        if seed not in seed_counts or (seed, row["row_hash"]) in identities:
            raise HeadFitError("row-universe seed identity differs")
        identities.add((seed, row["row_hash"]))
        seed_counts[seed] += 1
        current = {
            "element_id": row["element_id"],
            "source_locus_group_id": row["source_locus_group_id"],
            "long_range_block_id": row["long_range_block_id"],
            "outer_fold": row["outer_fold"],
            "unit_hash": row["unit_hash"],
            "block_hash": row["block_hash"],
        }
        previous = metadata.setdefault(row["element_id"], current)
        if previous != current:
            raise HeadFitError("element metadata differs across contexts or seeds")
    if set(seed_counts.values()) != {2066} or len(metadata) != 1033:
        raise HeadFitError("row-universe seed or element census differs")
    elements = sorted(metadata)
    if len({metadata[element]["long_range_block_id"] for element in elements}) != 239:
        raise HeadFitError("long-range block census differs")
    return rows, elements, metadata


def raw_binding(config: Mapping[str, Any], model_id: str) -> Mapping[str, Any]:
    matches = [raw for raw in config["raw_authority"] if model_id in raw["model_ids"]]
    if len(matches) != 1:
        raise HeadFitError(f"raw authority differs for {model_id}")
    return matches[0]


def load_raw_embeddings(
    root: Path,
    config: Mapping[str, Any],
    model_id: str,
    elements: Sequence[str],
) -> tuple[np.ndarray, list[dict[str, str]]]:
    binding = raw_binding(config, model_id)
    tree = safe_path(root, binding["tree_path"], label=f"{model_id} raw")
    fixture_tree = safe_path(
        root, binding["fixture_tree_path"], label=f"{model_id} fixture"
    )
    if file_sha256(tree / "ARTIFACTS.json") != binding["artifacts_sha256"]:
        raise HeadFitError(f"{model_id} raw ARTIFACTS hash differs")
    if file_sha256(fixture_tree / "ARTIFACTS.json") != binding["fixture_artifacts_sha256"]:
        raise HeadFitError(f"{model_id} fixture ARTIFACTS hash differs")
    verify_frozen_tree(tree)
    verify_frozen_tree(fixture_tree)
    fixture_rows = read_tsv(fixture_tree / binding["fixture_manifest_member"])
    fixture_by_id = {row["fixture_id"]: row for row in fixture_rows}
    with np.load(tree / f"raw/{model_id}/allele_embeddings.npz", allow_pickle=False) as data:
        fixture_ids = data["fixture_ids"].astype(str)
        allele_order = tuple(data["allele_order"].astype(str).tolist())
        embeddings = data["embeddings"]
    if (
        fixture_ids.shape != (1033,)
        or len(fixture_by_id) != 1033
        or set(fixture_ids) != set(fixture_by_id)
        or allele_order != ("REF", "ALT", "REF_RC", "ALT_RC")
        or embeddings.shape[:2] != (1033, 4)
        or embeddings.shape[2] < PROJECTION_WIDTH
        or not np.isfinite(embeddings).all()
    ):
        raise HeadFitError(f"{model_id} raw embedding contract differs")
    element_to_index = {
        fixture_by_id[fixture_id]["element_id"]: index
        for index, fixture_id in enumerate(fixture_ids)
    }
    if set(element_to_index) != set(elements):
        raise HeadFitError(f"{model_id} raw element universe differs")
    order = np.asarray([element_to_index[element] for element in elements], dtype=np.int64)
    fixture_by_element = {
        row["element_id"]: row for row in fixture_rows
    }
    return embeddings[order], [fixture_by_element[element] for element in elements]


def fit_inner_projection(
    embeddings: np.ndarray, training_mask: np.ndarray
) -> dict[str, np.ndarray]:
    values = np.asarray(embeddings, dtype=np.float64)
    mask = np.asarray(training_mask, dtype=bool)
    if (
        values.ndim != 3
        or values.shape[1] != 4
        or mask.shape != (values.shape[0],)
        or int(mask.sum()) < 2
        or not np.isfinite(values).all()
    ):
        raise HeadFitError("inner projection request differs")
    training = values[mask].reshape(-1, values.shape[2])
    mean = training.mean(axis=0)
    standard_deviation = np.maximum(training.std(axis=0), STD_FLOOR)
    normalized = (training - mean) / standard_deviation
    _, singular_values, right_vectors = np.linalg.svd(normalized, full_matrices=False)
    components = right_vectors[:PROJECTION_WIDTH].copy()
    for index in range(PROJECTION_WIDTH):
        pivot = int(np.argmax(np.abs(components[index])))
        if components[index, pivot] < 0:
            components[index] *= -1
    eigenvalues = singular_values[:PROJECTION_WIDTH] ** 2 / max(1, training.shape[0] - 1)
    eigenvalue_floor = max(
        EIGENVALUE_ABSOLUTE_FLOOR,
        EIGENVALUE_RELATIVE_FLOOR * float(eigenvalues[0]),
    )
    whitening_scale = np.sqrt(np.maximum(eigenvalues, eigenvalue_floor))
    return {
        "mean": mean,
        "standard_deviation": standard_deviation,
        "components": components,
        "whitening_scale": whitening_scale,
        "training_elements": np.asarray(mask.sum(), dtype=np.int64),
    }


def load_outer_features(
    projection_root: Path,
    model_id: str,
    held_fold: int,
    elements: Sequence[str],
    folds: np.ndarray,
) -> np.ndarray:
    path = projection_root / (
        f"open_sequence_projection/{model_id}/heldout_fold{held_fold}/head_features.npz"
    )
    with np.load(path, allow_pickle=False) as data:
        saved_elements = data["element_ids"].astype(str)
        saved_folds = data["outer_folds"].astype(np.int64)
        block_order = tuple(data["feature_block_order"].astype(str).tolist())
        features = data["features"]
    index = {element: position for position, element in enumerate(saved_elements)}
    if (
        len(index) != 1033
        or set(index) != set(elements)
        or block_order != ("REF", "ALT", "ALT_minus_REF", "absolute_ALT_minus_REF")
        or features.shape != (1033, 1024)
        or not np.isfinite(features).all()
    ):
        raise HeadFitError(f"{model_id} outer feature contract differs")
    order = np.asarray([index[element] for element in elements], dtype=np.int64)
    if not np.array_equal(saved_folds[order], folds):
        raise HeadFitError(f"{model_id} outer feature fold map differs")
    return features[order].astype(np.float64)


def standardization(values: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = values[mask].mean(axis=0)
    scale = np.maximum(values[mask].std(axis=0), 1.0e-6)
    return mean, scale


def fit_ridge_outer(
    *,
    inner_features: np.ndarray,
    outer_features: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    held_fold: int,
    alphas: Sequence[float],
) -> tuple[np.ndarray, dict[str, Any], dict[str, np.ndarray]]:
    inner_fold = (held_fold + 1) % 5
    outer_train = folds != held_fold
    outer_test = folds == held_fold
    inner_train = outer_train & (folds != inner_fold)
    inner_valid = folds == inner_fold
    inner_mean, inner_scale = standardization(inner_features, inner_train)
    x_train = (inner_features[inner_train] - inner_mean) / inner_scale
    x_valid = (inner_features[inner_valid] - inner_mean) / inner_scale
    scores: list[tuple[float, float]] = []
    for alpha in alphas:
        model = Ridge(alpha=float(alpha), fit_intercept=True, solver="lsqr", tol=1.0e-7)
        model.fit(x_train, outcomes[inner_train])
        prediction = model.predict(x_valid)
        scores.append((float(alpha), float(np.sqrt(np.mean((outcomes[inner_valid] - prediction) ** 2)))))
    constant_rmse = float(
        np.sqrt(np.mean((outcomes[inner_valid] - outcomes[inner_train].mean()) ** 2))
    )
    selected_alpha, selected_rmse = min(
        [*scores, (math.inf, constant_rmse)], key=lambda value: (value[1], value[0])
    )
    outer_mean, outer_scale = standardization(outer_features, outer_train)
    if math.isinf(selected_alpha):
        coefficient = np.zeros(outer_features.shape[1], dtype=np.float64)
        intercept = float(outcomes[outer_train].mean())
        prediction = np.full(int(outer_test.sum()), intercept, dtype=np.float64)
    else:
        final = Ridge(
            alpha=selected_alpha, fit_intercept=True, solver="lsqr", tol=1.0e-7
        )
        final.fit(
            (outer_features[outer_train] - outer_mean) / outer_scale,
            outcomes[outer_train],
        )
        prediction = final.predict(
            (outer_features[outer_test] - outer_mean) / outer_scale
        )
        coefficient = np.asarray(final.coef_, dtype=np.float64)
        intercept = float(final.intercept_)
    receipt = {
        "held_out_fold": held_fold,
        "inner_validation_fold": inner_fold,
        "inner_training_elements": int(inner_train.sum()),
        "inner_validation_elements": int(inner_valid.sum()),
        "outer_training_elements": int(outer_train.sum()),
        "outer_test_elements": int(outer_test.sum()),
        "selected_alpha": "constant_training_mean" if math.isinf(selected_alpha) else selected_alpha,
        "selected_inner_rmse": selected_rmse,
        "inner_rmse_by_alpha": {
            **{format(alpha, ".17g"): rmse for alpha, rmse in scores},
            "constant_training_mean": constant_rmse,
        },
        "held_out_outcomes_read_during_fit": False,
    }
    state = {
        "feature_mean": outer_mean,
        "feature_scale": outer_scale,
        "coefficient": coefficient,
        "intercept": np.asarray(intercept, dtype=np.float64),
        "selected_alpha": np.asarray(selected_alpha, dtype=np.float64),
    }
    return np.asarray(prediction, dtype=np.float64), receipt, state


def _mlp_class(input_width: int, hidden_width: int) -> Any:
    import torch
    from torch import nn

    class Head(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.layer1 = nn.Linear(input_width, hidden_width)
            self.activation = nn.GELU()
            self.layer2 = nn.Linear(hidden_width, 1)

        def forward(self, value: Any) -> Any:
            return self.layer2(self.activation(self.layer1(value))).squeeze(1)

    return Head


def select_mlp_epochs(
    *,
    features: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    held_fold: int,
    recipe: Mapping[str, Any],
) -> tuple[int, dict[str, Any]]:
    import torch

    inner_fold = (held_fold + 1) % 5
    inner_train = (folds != held_fold) & (folds != inner_fold)
    inner_valid = folds == inner_fold
    mean, scale = standardization(features, inner_train)
    outcome_mean = float(outcomes[inner_train].mean())
    outcome_scale = max(float(outcomes[inner_train].std()), 1.0e-6)
    x_train = torch.from_numpy(((features[inner_train] - mean) / scale).astype(np.float32))
    y_train = torch.from_numpy(((outcomes[inner_train] - outcome_mean) / outcome_scale).astype(np.float32))
    x_valid = torch.from_numpy(((features[inner_valid] - mean) / scale).astype(np.float32))
    y_valid = torch.from_numpy(((outcomes[inner_valid] - outcome_mean) / outcome_scale).astype(np.float32))
    seed = int(recipe["epoch_selection_seed"])
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    Head = _mlp_class(features.shape[1], int(recipe["hidden_width"]))
    model = Head()
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(recipe["learning_rate"]),
        weight_decay=float(recipe["weight_decay"]),
    )
    best_epoch, best_loss, stale = 1, math.inf, 0
    for epoch in range(1, int(recipe["maximum_epochs"]) + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        loss = torch.mean((model(x_train) - y_train) ** 2)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), float(recipe["gradient_clip_norm"]))
        optimizer.step()
        model.eval()
        with torch.no_grad():
            valid_loss = float(torch.mean((model(x_valid) - y_valid) ** 2).item())
        if valid_loss < best_loss - float(recipe["minimum_delta"]):
            best_epoch, best_loss, stale = epoch, valid_loss, 0
        else:
            stale += 1
        if stale >= int(recipe["patience"]):
            break
    receipt = {
        "held_out_fold": held_fold,
        "inner_validation_fold": inner_fold,
        "selection_seed": seed,
        "selected_epochs": best_epoch,
        "best_inner_standardized_mse": best_loss,
        "epochs_executed": epoch,
        "inner_training_elements": int(inner_train.sum()),
        "inner_validation_elements": int(inner_valid.sum()),
        "held_out_outcomes_read_during_fit": False,
    }
    return best_epoch, receipt


def fit_mlp_outer(
    *,
    features: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    held_fold: int,
    selected_epochs: int,
    seed: int,
    recipe: Mapping[str, Any],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    import torch

    train = folds != held_fold
    test = folds == held_fold
    mean, scale = standardization(features, train)
    outcome_mean = float(outcomes[train].mean())
    outcome_scale = max(float(outcomes[train].std()), 1.0e-6)
    x_train = torch.from_numpy(((features[train] - mean) / scale).astype(np.float32))
    y_train = torch.from_numpy(((outcomes[train] - outcome_mean) / outcome_scale).astype(np.float32))
    x_test = torch.from_numpy(((features[test] - mean) / scale).astype(np.float32))
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    Head = _mlp_class(features.shape[1], int(recipe["hidden_width"]))
    model = Head()
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(recipe["learning_rate"]),
        weight_decay=float(recipe["weight_decay"]),
    )
    for _ in range(selected_epochs):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        loss = torch.mean((model(x_train) - y_train) ** 2)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), float(recipe["gradient_clip_norm"]))
        optimizer.step()
    model.eval()
    with torch.no_grad():
        prediction = model(x_test).numpy().astype(np.float64) * outcome_scale + outcome_mean
    state: dict[str, np.ndarray] = {
        "feature_mean": mean,
        "feature_scale": scale,
        "outcome_mean": np.asarray(outcome_mean, dtype=np.float64),
        "outcome_scale": np.asarray(outcome_scale, dtype=np.float64),
        "selected_epochs": np.asarray(selected_epochs, dtype=np.int64),
        "seed": np.asarray(seed, dtype=np.int64),
    }
    for name, value in model.state_dict().items():
        state[name.replace(".", "__")] = value.detach().cpu().numpy()
    return prediction, state


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)


def fit_campaign(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config_path = arguments.config.resolve(strict=True)
    config = load_config(config_path)
    try:
        validate_config(config)
    except TaskSpecError as error:
        raise HeadFitError(str(error)) from error
    verify_taskspec(
        arguments.task_spec_tree.resolve(strict=True),
        arguments.task_spec_sha256,
        file_sha256(config_path),
    )
    if arguments.output.exists():
        raise HeadFitError("head-fit output exists")
    row_root = safe_path(root, config["row_authority"]["tree_path"], label="row authority")
    outcome_root = safe_path(
        root, config["outcome_authority"]["tree_path"], label="outcome authority"
    )
    projection_root = safe_path(
        root, config["projection_authority"]["tree_path"], label="projection authority"
    )
    for tree, authority, label in (
        (row_root, config["row_authority"], "row"),
        (outcome_root, config["outcome_authority"], "outcomes"),
        (projection_root, config["projection_authority"], "projection"),
    ):
        if file_sha256(tree / "ARTIFACTS.json") != authority["artifacts_sha256"]:
            raise HeadFitError(f"{label} ARTIFACTS hash differs")
        verify_frozen_tree(tree)
    row_path = row_root / config["row_authority"]["row_member"]
    outcome_path = outcome_root / config["outcome_authority"]["outcome_member"]
    rows, elements, metadata = load_row_universe(row_path)
    element_index = {element: index for index, element in enumerate(elements)}
    folds = np.asarray(
        [int(metadata[element]["outer_fold"].removeprefix("fold-")) for element in elements],
        dtype=np.int64,
    )
    if set(folds.tolist()) != set(range(5)):
        raise HeadFitError("outer fold coverage differs")
    outcomes = load_outcomes(outcome_path, set(elements))
    target = {
        context: np.asarray([outcomes[(element, context)]["mean"] for element in elements])
        for context in CONTEXTS
    }
    if any(not np.isfinite(values).all() for values in target.values()):
        raise HeadFitError("aggregated outcome is non-finite")
    raw: dict[str, np.ndarray] = {}
    fixture: dict[str, list[dict[str, str]]] = {}
    for model in MODELS:
        raw[model], fixture[model] = load_raw_embeddings(root, config, model, elements)
    allele_features = np.zeros((len(elements), 16), dtype=np.float64)
    bases = {base: index for index, base in enumerate("ACGT")}
    for index, row in enumerate(fixture["dnabert2"]):
        ref, alt = row["ref"], row["alt"]
        if ref not in bases or alt not in bases or ref == alt:
            raise HeadFitError("allele identity differs")
        allele_features[index, 4 * bases[ref] + bases[alt]] = 1.0

    arguments.output.mkdir(parents=True, mode=0o750)
    (arguments.output / "heads").mkdir()
    (arguments.output / "selection_projections").mkdir()
    outer_features: dict[tuple[str, int], np.ndarray] = {}
    inner_features: dict[tuple[str, int], np.ndarray] = {}
    projection_receipts = []
    for model in MODELS:
        model_projection_root = arguments.output / f"selection_projections/{model}"
        model_projection_root.mkdir()
        for held_fold in range(5):
            inner_fold = (held_fold + 1) % 5
            inner_train = (folds != held_fold) & (folds != inner_fold)
            parameters = fit_inner_projection(raw[model], inner_train)
            inner_feature = head_features(apply_projection(raw[model], parameters)).astype(
                np.float64
            )
            outer_feature = load_outer_features(
                projection_root, model, held_fold, elements, folds
            )
            inner_features[(model, held_fold)] = inner_feature
            outer_features[(model, held_fold)] = outer_feature
            np.savez_compressed(
                model_projection_root / f"outer{held_fold}_inner_projection.npz",
                **parameters,
            )
            projection_receipts.append(
                {
                    "model_id": model,
                    "outer_fold": held_fold,
                    "inner_validation_fold": inner_fold,
                    "inner_training_elements": int(inner_train.sum()),
                    "inner_projection_features_sha256": array_sha256(inner_feature),
                    "outer_projection_features_sha256": array_sha256(outer_feature),
                    "held_or_inner_validation_outcomes_used": False,
                }
            )

    predictions: dict[tuple[str, str, int, str], np.ndarray] = {}
    selection_rows: list[dict[str, Any]] = []
    alphas = tuple(float(value) for value in config["ridge_recipe"]["alpha_grid"])
    recipe = config["two_layer_recipe"]
    for context in CONTEXTS:
        y = target[context]
        zero = np.zeros(len(elements), dtype=np.float64)
        training_mean = np.empty(len(elements), dtype=np.float64)
        allele_ridge = np.empty(len(elements), dtype=np.float64)
        for held_fold in range(5):
            test = folds == held_fold
            train = ~test
            training_mean[test] = float(y[train].mean())
            ridge_prediction, selection, state = fit_ridge_outer(
                inner_features=allele_features,
                outer_features=allele_features,
                outcomes=y,
                folds=folds,
                held_fold=held_fold,
                alphas=alphas,
            )
            allele_ridge[test] = ridge_prediction
            state_root = arguments.output / f"heads/available_simple_controls/allele_identity_ridge/{context}"
            state_root.mkdir(parents=True, exist_ok=True)
            for seed in SEEDS:
                np.savez_compressed(state_root / f"seed{seed}_fold{held_fold}.npz", **state)
            selection_rows.append(
                {
                    "model_id": "available_simple_controls",
                    "head_id": "allele_identity_ridge",
                    "assay_context_id": context,
                    "seed": "deterministic_all_five",
                    **selection,
                }
            )
        for seed in SEEDS:
            predictions[("available_simple_controls", "zero", seed, context)] = zero.copy()
            predictions[("available_simple_controls", "training_mean", seed, context)] = training_mean.copy()
            predictions[("available_simple_controls", "allele_identity_ridge", seed, context)] = allele_ridge.copy()

        for model in MODELS:
            ridge_outputs = {
                "delta_ridge": np.empty(len(elements), dtype=np.float64),
                "full_ridge": np.empty(len(elements), dtype=np.float64),
            }
            mlp_outputs = {
                seed: np.empty(len(elements), dtype=np.float64) for seed in SEEDS
            }
            for held_fold in range(5):
                test = folds == held_fold
                inner_all = inner_features[(model, held_fold)]
                outer_all = outer_features[(model, held_fold)]
                slices = {
                    "delta_ridge": slice(512, 768),
                    "full_ridge": slice(0, 1024),
                }
                for head_id, feature_slice in slices.items():
                    ridge_prediction, selection, state = fit_ridge_outer(
                        inner_features=inner_all[:, feature_slice],
                        outer_features=outer_all[:, feature_slice],
                        outcomes=y,
                        folds=folds,
                        held_fold=held_fold,
                        alphas=alphas,
                    )
                    ridge_outputs[head_id][test] = ridge_prediction
                    state_root = arguments.output / f"heads/{model}/{head_id}/{context}"
                    state_root.mkdir(parents=True, exist_ok=True)
                    for seed in SEEDS:
                        np.savez_compressed(
                            state_root / f"seed{seed}_fold{held_fold}.npz", **state
                        )
                    selection_rows.append(
                        {
                            "model_id": model,
                            "head_id": head_id,
                            "assay_context_id": context,
                            "seed": "deterministic_all_five",
                            **selection,
                        }
                    )
                selected_epochs, selection = select_mlp_epochs(
                    features=inner_all,
                    outcomes=y,
                    folds=folds,
                    held_fold=held_fold,
                    recipe=recipe,
                )
                selection_rows.append(
                    {
                        "model_id": model,
                        "head_id": "two_layer_gelu",
                        "assay_context_id": context,
                        "seed": "selection_seed_1103_reused",
                        **selection,
                    }
                )
                for seed in SEEDS:
                    prediction, state = fit_mlp_outer(
                        features=outer_all,
                        outcomes=y,
                        folds=folds,
                        held_fold=held_fold,
                        selected_epochs=selected_epochs,
                        seed=seed,
                        recipe=recipe,
                    )
                    mlp_outputs[seed][test] = prediction
                    state_root = arguments.output / f"heads/{model}/two_layer_gelu/{context}"
                    state_root.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(
                        state_root / f"seed{seed}_fold{held_fold}.npz", **state
                    )
            for head_id, values in ridge_outputs.items():
                if not np.isfinite(values).all():
                    raise HeadFitError("ridge OOF prediction is incomplete")
                for seed in SEEDS:
                    predictions[(model, head_id, seed, context)] = values.copy()
            for seed, values in mlp_outputs.items():
                if not np.isfinite(values).all():
                    raise HeadFitError("two-layer OOF prediction is incomplete")
                predictions[(model, "two_layer_gelu", seed, context)] = values

    candidate_roster = [
        ("available_simple_controls", head) for head in BASELINES
    ] + [(model, head) for model in MODELS for head in HEADS]
    prediction_rows = []
    for model_id, head_id in candidate_roster:
        for row in rows:
            seed = int(row["seed"])
            context = row["assay_context_id"]
            index = element_index[row["element_id"]]
            prediction = float(predictions[(model_id, head_id, seed, context)][index])
            observed = float(target[context][index])
            if not math.isfinite(prediction) or not math.isfinite(observed):
                raise HeadFitError("standardized prediction row is non-finite")
            prediction_rows.append(
                {
                    "seed": seed,
                    "row_hash": row["row_hash"],
                    "unit_hash": row["unit_hash"],
                    "block_hash": row["block_hash"],
                    "stratum": row["stratum"],
                    "outer_fold": row["outer_fold"],
                    "study_id": row["study_id"],
                    "assay_context_id": context,
                    "element_id": row["element_id"],
                    "source_locus_group_id": row["source_locus_group_id"],
                    "long_range_block_id": row["long_range_block_id"],
                    "model_id": model_id,
                    "head_id": head_id,
                    "observed": format(observed, ".17g"),
                    "prediction": format(prediction, ".17g"),
                    "experimental_replicates": 4,
                    "biological_donors": 0,
                    "outcome_role": "exposed_development_MPRA_only",
                }
            )
    expected_rows = 12 * 10330
    if len(prediction_rows) != expected_rows:
        raise HeadFitError("standardized prediction denominator differs")
    write_gzip_tsv(
        arguments.output / "oof_predictions.tsv.gz", PREDICTION_FIELDS, prediction_rows
    )
    selection_fields = []
    for row in selection_rows:
        for field in row:
            if field not in selection_fields:
                selection_fields.append(field)
    write_tsv(arguments.output / "head_selection.tsv", selection_fields, selection_rows)
    projection_fields = tuple(projection_receipts[0])
    write_tsv(
        arguments.output / "inner_projection_audit.tsv",
        projection_fields,
        projection_receipts,
    )
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_prespecified_exposed_development_oof_fit",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "candidate_count": 12,
        "prediction_rows": expected_rows,
        "rows_per_candidate": 10330,
        "selected_elements": 1033,
        "contexts": list(CONTEXTS),
        "fixed_seeds": list(SEEDS),
        "outer_folds": 5,
        "long_range_blocks": 239,
        "inner_projection_refits": len(projection_receipts),
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "outcome_aggregation": config["task_spec"]["replicate_aggregation"],
        "outcome_role": "exposed_development_MPRA_only",
        "campaign_scope": "partial_common_head_campaign",
        "mandatory_baselines_complete": False,
        "mandatory_baselines_missing": list(MISSING_BASELINES),
        "shortlist_blocked": True,
        "finalist_claim_blocked": True,
        "complementarity_blocked": True,
        "conditional_trigger_blocked": True,
        "task_spec_artifacts_sha256": arguments.task_spec_sha256,
        "projection_artifacts_sha256": config["projection_authority"]["artifacts_sha256"],
        "outcome_artifacts_sha256": config["outcome_authority"]["artifacts_sha256"],
        "metrics_calculated": False,
        "sealed_assets_read": False,
        "champion_claim": False,
        "stack_fit": False,
        "residual_correlation_calculated": False,
        "conditional_model_built_or_fit": False,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-spec-tree", type=Path, required=True)
    parser.add_argument("--task-spec-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    fit_campaign(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
