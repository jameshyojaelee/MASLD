#!/usr/bin/env python3
"""Fit frozen common heads as resumable study-held-out shards and aggregate them.

The outcome-blind embedding output file supplies model, representation, checkpoint,
and exposure identity. A shard fits one head/seed/outer-fold combination and
emits predictions without observed labels. Aggregation never reads labels and
requires the complete frozen shard census.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import random
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
HEADS = ("linear", "two_layer_mlp")
SCREEN_SEEDS = (1103, 1201, 1301)
OUTER_FOLDS = tuple(range(5))
EXPOSURE_STATUSES = frozenset(
    {
        "clean_declared",
        "target_label_unexposed",
        "encoder_seen",
        "continual_seen",
        "reference_only",
        "downstream_demo",
        "unknown",
    }
)
BATCH_SIZE = 1024
FIXED_EPOCHS = 30
LEARNING_RATE = 0.001
WEIGHT_DECAY = 0.01
PREDICTION_FIELDS = (
    "row_id",
    "donor_id",
    "dataset",
    "outer_fold",
    "predicted_class",
    *(f"probability::{label}" for label in ROSTER),
)
FORBIDDEN_PREDICTION_FIELDS = frozenset(
    {
        "broad_label",
        "observed_label",
        "truth",
        "fibrosis",
        "fibrosis_score",
        "nas",
        "nas_score",
        "histology",
        "sealed_outcome",
    }
)
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_REPRESENTATION_FILE = re.compile(r"^(?P<representation>.+)_embeddings\.npz$")


class CommonHeadError(ValueError):
    """Raised when inputs, training, or predictions violate the common lane."""


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CommonHeadError(f"TSV lacks a header: {path}")
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


def _require_sha256(value: Any, label: str) -> str:
    text = str(value)
    if _SHA256.fullmatch(text) is None:
        raise CommonHeadError(f"{label} is not a lowercase SHA-256")
    return text


def _resolve_inside(root: Path, relative: Path, *, label: str) -> Path:
    if relative.is_absolute() or relative == Path("."):
        raise CommonHeadError(f"{label} must be a nonempty relative path")
    try:
        resolved = (root / relative).resolve(strict=True)
        resolved.relative_to(root.resolve(strict=True))
    except (OSError, ValueError) as error:
        raise CommonHeadError(f"{label} escapes its frozen root") from error
    return resolved


def derive_embedding_identity(
    embeddings_root: Path, embeddings_relative: Path
) -> tuple[dict[str, Any], Path]:
    """Derive all scientific identity from a verified frozen embedding output file."""

    manifest = verify_frozen_tree(embeddings_root)
    metadata = manifest.get("metadata", {})
    embedding_path = _resolve_inside(
        embeddings_root, embeddings_relative, label="embedding path"
    )
    match = _REPRESENTATION_FILE.fullmatch(embedding_path.name)
    if match is None:
        raise CommonHeadError("embedding filename does not bind a representation")
    representation_id = match.group("representation")
    receipt_path = embedding_path.parent / "receipt.json"
    if not receipt_path.is_file():
        raise CommonHeadError("frozen embedding receipt is absent")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CommonHeadError("frozen embedding receipt is invalid") from error

    model_id = metadata.get("model_id")
    checkpoint_sha256 = _require_sha256(
        receipt.get("checkpoint_sha256"), "embedding checkpoint"
    )
    exposure_status = metadata.get("exposure_status")
    policies = receipt.get("policies")
    if (
        not isinstance(model_id, str)
        or not model_id
        or receipt.get("model_id", model_id) != model_id
        or exposure_status not in EXPOSURE_STATUSES
        or receipt.get("exposure_status") != exposure_status
        or not isinstance(policies, (dict, list))
        or representation_id not in policies
        or metadata.get("dataset_view_id") != DATASET_VIEW_ID
        or metadata.get("development_rows") != 50_000
        or metadata.get("evaluation_labels_read") is not False
        or metadata.get("sealed_outcomes_read") is not False
        or receipt.get("rows") != 50_000
        or receipt.get("evaluation_label_columns_read") not in ([], False)
        or receipt.get("sealed_outcomes_read") is not False
    ):
        raise CommonHeadError("frozen embedding identity or firewall differs")
    width = receipt.get("embedding_width")
    if not isinstance(width, int) or width < 1:
        raise CommonHeadError("embedding width is not frozen")
    identity = {
        "model_id": model_id,
        "representation_id": representation_id,
        "checkpoint_sha256": checkpoint_sha256,
        "exposure_status": exposure_status,
        "embedding_width": width,
        "embedding_relative": embeddings_relative.as_posix(),
        "sealed_champion_eligible": bool(
            metadata.get("sealed_champion_eligible", False)
        ),
    }
    if identity["sealed_champion_eligible"] and exposure_status not in {
        "clean_declared",
        "target_label_unexposed",
    }:
        raise CommonHeadError("embedding champion eligibility contradicts exposure")
    return identity, embedding_path


def donor_class_weights(
    donors: Sequence[str], labels: Sequence[int], *, require_full_roster: bool = True
) -> np.ndarray:
    if len(donors) != len(labels) or not len(donors):
        raise CommonHeadError("donor and label vectors differ")
    counts: dict[tuple[str, int], int] = {}
    donors_by_class: dict[int, set[str]] = {}
    for donor, label in zip(donors, labels, strict=True):
        key = (str(donor), int(label))
        counts[key] = counts.get(key, 0) + 1
        donors_by_class.setdefault(int(label), set()).add(str(donor))
    if require_full_roster and set(donors_by_class) != set(range(len(ROSTER))):
        raise CommonHeadError("training partition lacks a frozen class")
    class_count = len(donors_by_class)
    weights = np.asarray(
        [
            1.0
            / (
                class_count
                * len(donors_by_class[int(label)])
                * counts[(str(donor), int(label))]
            )
            for donor, label in zip(donors, labels, strict=True)
        ],
        dtype=np.float32,
    )
    return weights / weights.sum() * len(weights)


def validate_study_outer_split(
    donors: np.ndarray, datasets: np.ndarray, outer: np.ndarray
) -> None:
    if not (donors.ndim == datasets.ndim == outer.ndim == 1) or not (
        len(donors) == len(datasets) == len(outer)
    ):
        raise CommonHeadError("study outer split arrays differ")
    if set(map(int, outer)) != set(OUTER_FOLDS):
        raise CommonHeadError("outer fold roster differs")
    if any(len(set(map(int, outer[donors == donor]))) != 1 for donor in set(donors)):
        raise CommonHeadError("one donor crosses outer folds")
    if any(
        len(set(map(int, outer[datasets == dataset]))) != 1
        for dataset in set(datasets)
    ):
        raise CommonHeadError("one study crosses outer folds")
    for held_outer in OUTER_FOLDS:
        training_studies = set(datasets[outer != held_outer])
        test_studies = set(datasets[outer == held_outer])
        if not test_studies or training_studies & test_studies:
            raise CommonHeadError("outer train/test studies are not disjoint")


def validate_inner_assignment(
    training_donors: np.ndarray, assignment: Mapping[str, int]
) -> np.ndarray:
    donor_roster = set(map(str, training_donors))
    if set(assignment) != donor_roster:
        raise CommonHeadError("inner assignment does not match outer training donors")
    if any(type(value) is not int or value not in OUTER_FOLDS for value in assignment.values()):
        raise CommonHeadError("inner fold is outside 0..4")
    row_folds = np.asarray([assignment[str(donor)] for donor in training_donors], dtype=np.int8)
    if set(map(int, row_folds)) != set(OUTER_FOLDS):
        raise CommonHeadError("an inner fold is empty")
    coverage = sum((row_folds == fold).astype(np.int8) for fold in OUTER_FOLDS)
    if not np.all(coverage == 1):
        raise CommonHeadError("inner rows do not have one-time OOF coverage")
    return row_folds


def build_head(head_id: str, width: int):
    from torch import nn

    if head_id == "linear":
        return nn.Linear(width, len(ROSTER))
    if head_id == "two_layer_mlp":
        hidden = min(256, width)
        return nn.Sequential(
            nn.Linear(width, hidden),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(hidden, len(ROSTER)),
        )
    raise CommonHeadError(f"unknown common head: {head_id}")


def _weighted_cross_entropy(logits, targets, weights):
    import torch

    losses = torch.nn.functional.cross_entropy(logits, targets, reduction="none")
    return (losses * weights).sum() / weights.sum()


def _fit_temperature(logits: np.ndarray, targets: np.ndarray, weights: np.ndarray) -> float:
    import torch

    values = torch.tensor(logits, dtype=torch.float64)
    outcomes = torch.tensor(targets, dtype=torch.long)
    sample_weights = torch.tensor(weights, dtype=torch.float64)
    log_temperature = torch.tensor(0.0, dtype=torch.float64, requires_grad=True)
    optimizer = torch.optim.LBFGS(
        [log_temperature], lr=0.1, max_iter=100, line_search_fn="strong_wolfe"
    )

    def closure():
        optimizer.zero_grad(set_to_none=True)
        temperature = log_temperature.exp().clamp(0.05, 20.0)
        loss = _weighted_cross_entropy(values / temperature, outcomes, sample_weights)
        loss.backward()
        return loss

    optimizer.step(closure)
    temperature = float(log_temperature.exp().clamp(0.05, 20.0).detach())
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise CommonHeadError("temperature calibration failed")
    return temperature


def _fit_fixed_head(
    *,
    head_id: str,
    features: np.ndarray,
    targets: np.ndarray,
    donors: np.ndarray,
    fit_indices: np.ndarray,
    prediction_indices: np.ndarray,
    seed: int,
):
    """Fit exactly FIXED_EPOCHS; prediction-row labels are never read."""

    import torch

    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))))

    fit_x = features[fit_indices].astype(np.float32, copy=False)
    prediction_x = features[prediction_indices].astype(np.float32, copy=False)
    mean = fit_x.mean(axis=0, dtype=np.float64).astype(np.float32)
    scale = fit_x.std(axis=0, dtype=np.float64).astype(np.float32)
    scale[scale == 0.0] = 1.0
    fit_x = (fit_x - mean) / scale
    prediction_x = (prediction_x - mean) / scale
    fit_y = targets[fit_indices]
    fit_donors = donors[fit_indices]
    fit_weights = donor_class_weights(fit_donors, fit_y)

    device = torch.device("cuda")
    model = build_head(head_id, features.shape[1]).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )
    x = torch.from_numpy(fit_x).to(device)
    y = torch.from_numpy(fit_y.astype(np.int64, copy=False)).to(device)
    w = torch.from_numpy(fit_weights).to(device)
    prediction_x_gpu = torch.from_numpy(prediction_x).to(device)
    generator = torch.Generator().manual_seed(seed)
    for _epoch in range(FIXED_EPOCHS):
        model.train()
        order = torch.randperm(len(x), generator=generator)
        for start in range(0, len(x), BATCH_SIZE):
            batch = order[start : start + BATCH_SIZE].to(device)
            optimizer.zero_grad(set_to_none=True)
            loss = _weighted_cross_entropy(model(x[batch]), y[batch], w[batch])
            loss.backward()
            optimizer.step()
    model.eval()
    with torch.no_grad():
        final_training_loss = float(_weighted_cross_entropy(model(x), y, w).item())
        prediction_logits = model(prediction_x_gpu).float().cpu().numpy()
    state = {
        key: value.detach().cpu().clone() for key, value in model.state_dict().items()
    }
    return model, mean, scale, prediction_logits, state, {
        "fixed_epochs": FIXED_EPOCHS,
        "final_training_weighted_cross_entropy": final_training_loss,
        "fit_rows": len(fit_indices),
        "prediction_rows": len(prediction_indices),
        "fit_donors": len(set(map(str, fit_donors))),
        "prediction_donors": len(set(map(str, donors[prediction_indices]))),
        "validation_labels_used_for_checkpoint_selection": False,
    }


def _save_model_state(
    path: Path,
    *,
    state: Mapping[str, Any],
    mean: np.ndarray,
    scale: np.ndarray,
) -> None:
    arrays = {"feature_mean": mean, "feature_scale": scale}
    arrays.update(
        {
            f"state::{key}": value.detach().cpu().numpy()
            for key, value in state.items()
        }
    )
    with path.open("xb") as handle:
        np.savez_compressed(handle, **arrays)


def _validate_cuda_runtime() -> None:
    import torch

    if os.environ.get("CUBLAS_WORKSPACE_CONFIG") not in {":4096:8", ":16:8"}:
        raise CommonHeadError("deterministic CUBLAS workspace is not configured")
    if os.environ.get("PYTHONHASHSEED") != "20260824":
        raise CommonHeadError("PYTHONHASHSEED differs")
    if (
        not torch.cuda.is_available()
        or torch.cuda.device_count() != 1
        or "L40S" not in torch.cuda.get_device_name(0)
        or list(torch.cuda.get_device_capability(0)) != [8, 9]
    ):
        raise CommonHeadError("common heads require one L40S sm89 device")


def _load_inputs(
    *,
    embeddings_root: Path,
    embeddings_relative: Path,
    source: Path,
    split: Path,
    expected_embeddings_artifacts_sha256: str,
    expected_source_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
) -> dict[str, Any]:
    expected = {
        embeddings_root / "ARTIFACTS.json": expected_embeddings_artifacts_sha256,
        source / "ARTIFACTS.json": expected_source_artifacts_sha256,
        split / "ARTIFACTS.json": expected_split_artifacts_sha256,
    }
    for path, digest in expected.items():
        if sha256_file(path) != _require_sha256(digest, f"expected hash for {path}"):
            raise CommonHeadError(f"input ARTIFACTS SHA-256 differs: {path}")
    identity, embedding_path = derive_embedding_identity(
        embeddings_root, embeddings_relative
    )
    source_manifest = verify_frozen_tree(source)
    split_manifest = verify_frozen_tree(split)
    if (
        source_manifest["metadata"].get("subset_id") != DATASET_VIEW_ID
        or source_manifest["metadata"].get("row_count") != 50_000
        or split_manifest["metadata"].get("split_id") != SPLIT_ID
        or split_manifest["metadata"].get("target_labels_used_for_assignment") is not False
        or split_manifest["metadata"].get("sealed_outcomes_read") is not False
    ):
        raise CommonHeadError("frozen source or split metadata differs")
    with np.load(embedding_path, allow_pickle=False) as bundle:
        if set(bundle.files) != {"embeddings", "outer_folds", "row_ids"}:
            raise CommonHeadError("embedding bundle fields differ")
        features = np.asarray(bundle["embeddings"], dtype=np.float32)
        activation_outer = np.asarray(bundle["outer_folds"], dtype=np.int8)
        embedding_row_ids = np.asarray(bundle["row_ids"]).astype(str)
    if (
        features.shape != (50_000, identity["embedding_width"])
        or not np.all(np.isfinite(features))
        or len(set(embedding_row_ids)) != 50_000
        or activation_outer.shape != (50_000,)
    ):
        raise CommonHeadError("embedding matrix shape or values differ")

    selection = _read_tsv(source / "selection.tsv")
    split_rows = _read_tsv(split / "row_outer_folds.tsv")
    inner_rows = _read_tsv(split / "inner_donor_folds.tsv")
    row_ids = np.asarray([row["row_id"] for row in selection], dtype=str)
    donors = np.asarray([row["donor_id"] for row in selection], dtype=str)
    datasets = np.asarray([row["dataset"] for row in selection], dtype=str)
    labels = np.asarray([row["broad_label"] for row in selection], dtype=str)
    if (
        embedding_row_ids.tolist() != row_ids.tolist()
        or [row["row_id"] for row in split_rows] != row_ids.tolist()
        or [row["donor_id"] for row in split_rows] != donors.tolist()
        or [row["dataset"] for row in split_rows] != datasets.tolist()
        or set(labels) != set(ROSTER)
    ):
        raise CommonHeadError("embedding, source-label, and split rows differ")
    outer = np.asarray([int(row["outer_fold"]) for row in split_rows], dtype=np.int8)
    validate_study_outer_split(donors, datasets, outer)
    targets = np.asarray([ROSTER.index(label) for label in labels], dtype=np.int64)
    inner_by_outer: dict[int, dict[str, int]] = {fold: {} for fold in OUTER_FOLDS}
    for row in inner_rows:
        held = int(row["held_outer_fold"])
        inner = int(row["inner_fold"])
        donor = row["donor_id"]
        if held not in OUTER_FOLDS or donor in inner_by_outer[held]:
            raise CommonHeadError("inner donor or held fold differs")
        inner_by_outer[held][donor] = inner
    for held_outer in OUTER_FOLDS:
        validate_inner_assignment(donors[outer != held_outer], inner_by_outer[held_outer])
    return {
        "identity": identity,
        "features": features,
        "activation_outer": activation_outer,
        "row_ids": row_ids,
        "donors": donors,
        "datasets": datasets,
        "targets": targets,
        "outer": outer,
        "inner_by_outer": inner_by_outer,
    }


def run_shard(
    *,
    embeddings_root: Path,
    embeddings_relative: Path,
    source: Path,
    split: Path,
    output: Path,
    expected_embeddings_artifacts_sha256: str,
    expected_source_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
    head_id: str,
    screen_seed: int,
    held_outer: int,
) -> None:
    import torch

    if output.exists() or head_id not in HEADS or screen_seed not in SCREEN_SEEDS:
        raise CommonHeadError("output, head, or screen seed differs")
    if held_outer not in OUTER_FOLDS:
        raise CommonHeadError("outer fold differs")
    _validate_cuda_runtime()
    inputs = _load_inputs(
        embeddings_root=embeddings_root,
        embeddings_relative=embeddings_relative,
        source=source,
        split=split,
        expected_embeddings_artifacts_sha256=expected_embeddings_artifacts_sha256,
        expected_source_artifacts_sha256=expected_source_artifacts_sha256,
        expected_split_artifacts_sha256=expected_split_artifacts_sha256,
    )
    features = inputs["features"]
    targets = inputs["targets"]
    donors = inputs["donors"]
    datasets = inputs["datasets"]
    outer = inputs["outer"]
    training = np.flatnonzero(outer != held_outer)
    test = np.flatnonzero(outer == held_outer)
    if set(datasets[training]) & set(datasets[test]):
        raise CommonHeadError("shard train/test studies overlap")
    row_inner_folds = validate_inner_assignment(
        donors[training], inputs["inner_by_outer"][held_outer]
    )
    output.mkdir(mode=0o750)
    state_root = output / "states"
    state_root.mkdir()
    oof_logits = np.full((len(training), len(ROSTER)), np.nan, dtype=np.float64)
    oof_coverage = np.zeros(len(training), dtype=np.int8)
    inner_records = []
    for inner_fold in OUTER_FOLDS:
        validation_local = row_inner_folds == inner_fold
        validation_indices = training[validation_local]
        fit_indices = training[~validation_local]
        seed = (
            screen_seed
            + held_outer * 101
            + inner_fold * 17
            + HEADS.index(head_id)
        )
        model, mean, scale, logits, state, record = _fit_fixed_head(
            head_id=head_id,
            features=features,
            targets=targets,
            donors=donors,
            fit_indices=fit_indices,
            prediction_indices=validation_indices,
            seed=seed,
        )
        _save_model_state(
            state_root / f"inner{inner_fold}.npz",
            state=state,
            mean=mean,
            scale=scale,
        )
        oof_logits[validation_local] = logits
        oof_coverage[validation_local] += 1
        inner_records.append({"inner_fold": inner_fold, "seed": seed, **record})
        del model
        torch.cuda.empty_cache()
    if not np.all(oof_coverage == 1) or not np.all(np.isfinite(oof_logits)):
        raise CommonHeadError("OOF logits lack exact one-time coverage")
    temperature = _fit_temperature(
        oof_logits,
        targets[training],
        donor_class_weights(donors[training], targets[training]),
    )
    final_seed = screen_seed + held_outer * 101 + 1009 + HEADS.index(head_id)
    final_model, final_mean, final_scale, final_logits, final_state, final_record = (
        _fit_fixed_head(
            head_id=head_id,
            features=features,
            targets=targets,
            donors=donors,
            fit_indices=training,
            prediction_indices=test,
            seed=final_seed,
        )
    )
    _save_model_state(
        state_root / "final_outer_training_refit.npz",
        state=final_state,
        mean=final_mean,
        scale=final_scale,
    )
    final_model.eval()
    with torch.no_grad():
        probabilities = (
            torch.softmax(torch.from_numpy(final_logits).to("cuda") / temperature, dim=1)
            .double()
            .cpu()
            .numpy()
        )
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    if (
        not np.all(np.isfinite(probabilities))
        or np.any(probabilities < 0.0)
        or not np.allclose(
            probabilities.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12
        )
    ):
        raise CommonHeadError("common-head probabilities differ")
    predicted = [ROSTER[int(np.argmax(row))] for row in probabilities]
    _write_tsv(
        output / "predictions.tsv",
        PREDICTION_FIELDS,
        (
            {
                "row_id": inputs["row_ids"][index],
                "donor_id": donors[index],
                "dataset": datasets[index],
                "outer_fold": int(outer[index]),
                "predicted_class": predicted[position],
                **{
                    f"probability::{label}": format(
                        float(probabilities[position, class_index]), ".17g"
                    )
                    for class_index, label in enumerate(ROSTER)
                },
            }
            for position, index in enumerate(test)
        ),
    )
    receipt = {
        "schema_version": "masld-bench-common-cell-head-shard-v1",
        "status": "pass_development_prediction_shard",
        **inputs["identity"],
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "head_id": head_id,
        "screen_seed": screen_seed,
        "outer_fold": held_outer,
        "training_rows": len(training),
        "training_donors": len(set(donors[training])),
        "training_studies": sorted(set(datasets[training])),
        "test_rows": len(test),
        "test_donors": len(set(donors[test])),
        "test_studies": sorted(set(datasets[test])),
        "fixed_epochs": FIXED_EPOCHS,
        "checkpoint_selection_uses_validation_labels": False,
        "calibration": "fixed_epoch_donor_cross_fitted_temperature",
        "temperature": temperature,
        "inner_models": inner_records,
        "inner_models_used_for_outer_prediction": False,
        "oof_rows": len(training),
        "oof_exact_one_time_coverage": True,
        "final_fit_rows": int(final_record["fit_rows"]),
        "final_fit_donors": int(final_record["fit_donors"]),
        "final_fit_studies": len(set(datasets[training])),
        "final_fit_seed": final_seed,
        "final_fit_fixed_epochs": int(final_record["fixed_epochs"]),
        "outer_prediction_model": "fixed_epoch_final_outer_training_refit",
        "outer_prediction_uses_final_refit": True,
        "activation_fold_annotations_used_for_head_fitting": False,
        "activation_vs_study_fold_mismatch_rows": int(
            np.sum(inputs["activation_outer"] != outer)
        ),
        "metrics_calculated": False,
        "development_labels_read": ["broad_label"],
        "prediction_tables_contain_observed_labels": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "input_artifacts_sha256": {
            "embeddings": expected_embeddings_artifacts_sha256,
            "source": expected_source_artifacts_sha256,
            "split": expected_split_artifacts_sha256,
        },
    }
    write_json_exclusive(output / "prediction_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "cell_foundation_common_head_prediction_shard",
            **inputs["identity"],
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "head_id": head_id,
            "screen_seed": screen_seed,
            "outer_fold": held_outer,
            "metrics_calculated": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )


def _read_prediction_shard(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != PREDICTION_FIELDS:
            raise CommonHeadError("prediction shard fields differ")
        if FORBIDDEN_PREDICTION_FIELDS & set(reader.fieldnames or ()):
            raise CommonHeadError("prediction shard exposes forbidden fields")
        rows = list(reader)
    for row in rows:
        try:
            probabilities = np.asarray(
                [float(row[f"probability::{label}"]) for label in ROSTER],
                dtype=np.float64,
            )
        except (KeyError, ValueError) as error:
            raise CommonHeadError("prediction shard probabilities are invalid") from error
        if (
            row["predicted_class"] not in ROSTER
            or not np.all(np.isfinite(probabilities))
            or np.any(probabilities < 0.0)
            or not np.isclose(probabilities.sum(), 1.0, rtol=0.0, atol=1.0e-12)
            or row["predicted_class"] != ROSTER[int(np.argmax(probabilities))]
        ):
            raise CommonHeadError("prediction shard values differ")
    return rows


def aggregate_shards(
    *,
    shard_manifest: Path,
    source: Path,
    split: Path,
    output: Path,
    expected_source_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
) -> None:
    if output.exists():
        raise CommonHeadError("refusing to overwrite common-head aggregate")
    for root, expected in (
        (source, expected_source_artifacts_sha256),
        (split, expected_split_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != _require_sha256(
            expected, f"expected hash for {root}"
        ):
            raise CommonHeadError("aggregate input ARTIFACTS SHA-256 differs")
        verify_frozen_tree(root)
    source_metadata = verify_frozen_tree(source)["metadata"]
    split_metadata = verify_frozen_tree(split)["metadata"]
    if (
        source_metadata.get("subset_id") != DATASET_VIEW_ID
        or source_metadata.get("row_count") != 50_000
        or split_metadata.get("split_id") != SPLIT_ID
        or split_metadata.get("target_labels_used_for_assignment") is not False
        or split_metadata.get("sealed_outcomes_read") is not False
    ):
        raise CommonHeadError("aggregate source or split metadata differs")
    split_rows = _read_tsv(split / "row_outer_folds.tsv")
    expected_rows = [row["row_id"] for row in split_rows]
    expected_outer = {row["row_id"]: int(row["outer_fold"]) for row in split_rows}
    expected_donor = {row["row_id"]: row["donor_id"] for row in split_rows}
    expected_dataset = {row["row_id"]: row["dataset"] for row in split_rows}
    manifest_rows = _read_tsv(shard_manifest)
    manifest_fields = tuple(manifest_rows[0].keys()) if manifest_rows else ()
    if manifest_fields != (
        "shard_root",
        "artifacts_sha256",
    ):
        raise CommonHeadError("shard manifest fields differ")
    shards: dict[tuple[str, int, int], tuple[dict[str, Any], list[dict[str, str]]]] = {}
    identity: dict[str, Any] | None = None
    for row in manifest_rows:
        shard_root = Path(row["shard_root"])
        if not shard_root.is_absolute() or not shard_root.is_dir():
            raise CommonHeadError("shard root is not an absolute directory")
        if sha256_file(shard_root / "ARTIFACTS.json") != _require_sha256(
            row["artifacts_sha256"], "shard artifact"
        ):
            raise CommonHeadError("shard ARTIFACTS SHA-256 differs")
        verify_frozen_tree(shard_root)
        receipt = json.loads(
            (shard_root / "prediction_receipt.json").read_text(encoding="utf-8")
        )
        if (
            receipt.get("status") != "pass_development_prediction_shard"
            or receipt.get("dataset_view_id") != DATASET_VIEW_ID
            or receipt.get("split_id") != SPLIT_ID
            or receipt.get("fixed_epochs") != FIXED_EPOCHS
            or receipt.get("checkpoint_selection_uses_validation_labels") is not False
            or receipt.get("oof_exact_one_time_coverage") is not True
            or receipt.get("inner_models_used_for_outer_prediction") is not False
            or receipt.get("final_fit_rows") != receipt.get("training_rows")
            or receipt.get("final_fit_donors") != receipt.get("training_donors")
            or receipt.get("final_fit_studies") != len(receipt.get("training_studies", []))
            or receipt.get("final_fit_fixed_epochs") != FIXED_EPOCHS
            or receipt.get("outer_prediction_uses_final_refit") is not True
            or receipt.get("outer_prediction_model")
            != "fixed_epoch_final_outer_training_refit"
            or receipt.get("prediction_tables_contain_observed_labels") is not False
            or receipt.get("sealed_outcomes_read") is not False
            or receipt.get("input_artifacts_sha256", {}).get("source")
            != expected_source_artifacts_sha256
            or receipt.get("input_artifacts_sha256", {}).get("split")
            != expected_split_artifacts_sha256
        ):
            raise CommonHeadError("prediction shard receipt differs")
        observed_identity = {
            key: receipt[key]
            for key in (
                "model_id",
                "representation_id",
                "checkpoint_sha256",
                "exposure_status",
                "embedding_width",
                "embedding_relative",
                "sealed_champion_eligible",
            )
        }
        observed_identity["embeddings_artifacts_sha256"] = receipt[
            "input_artifacts_sha256"
        ]["embeddings"]
        if identity is None:
            identity = observed_identity
        elif observed_identity != identity:
            raise CommonHeadError("shard identities differ")
        key = (
            str(receipt["head_id"]),
            int(receipt["screen_seed"]),
            int(receipt["outer_fold"]),
        )
        if key in shards:
            raise CommonHeadError("duplicate common-head shard")
        predictions = _read_prediction_shard(shard_root / "predictions.tsv")
        if any(int(item["outer_fold"]) != key[2] for item in predictions):
            raise CommonHeadError("prediction rows differ from shard outer fold")
        shards[key] = (receipt, predictions)
    expected_keys = {
        (head, seed, fold)
        for head in HEADS
        for seed in SCREEN_SEEDS
        for fold in OUTER_FOLDS
    }
    if set(shards) != expected_keys or identity is None:
        raise CommonHeadError("common-head shard census is incomplete")

    output.mkdir(mode=0o750)
    prediction_root = output / "predictions"
    prediction_root.mkdir()
    row_position = {row_id: index for index, row_id in enumerate(expected_rows)}
    for head in HEADS:
        for seed in SCREEN_SEEDS:
            combined = [
                item
                for fold in OUTER_FOLDS
                for item in shards[(head, seed, fold)][1]
            ]
            if (
                len(combined) != len(expected_rows)
                or len({item["row_id"] for item in combined}) != len(expected_rows)
                or set(item["row_id"] for item in combined) != set(expected_rows)
                or any(
                    int(item["outer_fold"]) != expected_outer[item["row_id"]]
                    or item["donor_id"] != expected_donor[item["row_id"]]
                    or item["dataset"] != expected_dataset[item["row_id"]]
                    for item in combined
                )
            ):
                raise CommonHeadError("aggregate row coverage differs")
            combined.sort(key=lambda item: row_position[item["row_id"]])
            _write_tsv(
                prediction_root / f"{head}__seed{seed}.tsv",
                PREDICTION_FIELDS,
                combined,
            )
    receipt = {
        "schema_version": "masld-bench-common-cell-head-study-screen-v2",
        "status": "pass_development_predictions",
        **identity,
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": len(expected_rows),
        "heads": list(HEADS),
        "screen_seeds": list(SCREEN_SEEDS),
        "outer_folds": len(OUTER_FOLDS),
        "logical_shards": len(expected_keys),
        "fixed_epochs": FIXED_EPOCHS,
        "checkpoint_selection_uses_validation_labels": False,
        "inner_models_used_for_outer_prediction": False,
        "outer_prediction_uses_final_refit": True,
        "metrics_calculated": False,
        "prediction_tables_contain_observed_labels": False,
        "development_labels_read_by_aggregator": [],
        "histology_read": False,
        "sealed_outcomes_read": False,
        "shard_manifest_sha256": sha256_file(shard_manifest),
    }
    write_json_exclusive(output / "prediction_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "cell_foundation_common_head_predictions",
            **identity,
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "rows": len(expected_rows),
            "logical_shards": len(expected_keys),
            "metrics_calculated": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    actions = value.add_subparsers(dest="action", required=True)
    shard = actions.add_parser("shard")
    shard.add_argument("--embeddings-root", required=True, type=Path)
    shard.add_argument("--embeddings-relative", required=True, type=Path)
    shard.add_argument("--source", required=True, type=Path)
    shard.add_argument("--split", required=True, type=Path)
    shard.add_argument("--output", required=True, type=Path)
    shard.add_argument("--expected-embeddings-artifacts-sha256", required=True)
    shard.add_argument("--expected-source-artifacts-sha256", required=True)
    shard.add_argument("--expected-split-artifacts-sha256", required=True)
    shard.add_argument("--head-id", required=True, choices=HEADS)
    shard.add_argument("--screen-seed", required=True, type=int, choices=SCREEN_SEEDS)
    shard.add_argument("--outer-fold", required=True, type=int, choices=OUTER_FOLDS)
    aggregate = actions.add_parser("aggregate")
    aggregate.add_argument("--shard-manifest", required=True, type=Path)
    aggregate.add_argument("--source", required=True, type=Path)
    aggregate.add_argument("--split", required=True, type=Path)
    aggregate.add_argument("--output", required=True, type=Path)
    aggregate.add_argument("--expected-source-artifacts-sha256", required=True)
    aggregate.add_argument("--expected-split-artifacts-sha256", required=True)
    return value


def main() -> int:
    arguments = parser().parse_args()
    if arguments.action == "shard":
        run_shard(
            embeddings_root=arguments.embeddings_root,
            embeddings_relative=arguments.embeddings_relative,
            source=arguments.source,
            split=arguments.split,
            output=arguments.output,
            expected_embeddings_artifacts_sha256=arguments.expected_embeddings_artifacts_sha256,
            expected_source_artifacts_sha256=arguments.expected_source_artifacts_sha256,
            expected_split_artifacts_sha256=arguments.expected_split_artifacts_sha256,
            head_id=arguments.head_id,
            screen_seed=arguments.screen_seed,
            held_outer=arguments.outer_fold,
        )
    else:
        aggregate_shards(
            shard_manifest=arguments.shard_manifest,
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
