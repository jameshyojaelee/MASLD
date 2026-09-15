#!/usr/bin/env python3
"""Fit one outer-fold bundle of task-native sequence controls and predict only."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import random
import time
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np
from safetensors.torch import save_file
import torch
from torch.nn import functional as F

from masld_bench.artifacts import verify_frozen_tree
from scripts.gse281364_task_native_sequence_controls import (
    CONTEXTS,
    MODEL_IDS,
    allele_pair_activity,
    build_model,
    joint_orientation_augmentation,
)


SCHEMA = "masld-bench-gse281364-task-native-sequence-production-v1"
FOLDS = tuple(f"fold-{index}" for index in range(5))
SEEDS = (1103, 2909, 4721, 6673, 8111)
CHANNEL_INDEX = {base: index for index, base in enumerate(("A", "G", "C", "T"))}


class SequenceTrainingError(RuntimeError):
    """Raised when a fit, split, or prediction-only requirement differs."""


def file_sha256(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    value = tomllib.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SequenceTrainingError("campaign config must be a table")
    return value


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SequenceTrainingError("JSON authority must be an object")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def validate_campaign(config: Mapping[str, Any]) -> None:
    execution = config.get("execution", {})
    output = config.get("output_contract", {})
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status")
        != "authorized_prediction_only_development_campaign"
        or config.get("model_training_authorized") is not True
        or config.get("prediction_generation_authorized") is not True
        or config.get("adapter_metric_calculation_authorized") is not False
        or config.get("adapter_scoring_authorized") is not False
        or config.get("model_ranking_authorized") is not False
        or config.get("promotion_authorized") is not False
        or tuple(execution.get("models", ())) != MODEL_IDS
        or tuple(execution.get("outer_folds", ())) != FOLDS
        or tuple(execution.get("seeds", ())) != SEEDS
        or execution.get("fit_count") != 50
        or execution.get("held_test_outcomes_visible_to_training_adapter")
        is not False
        or output.get("checkpoint_format") != "safetensors_tensor_only"
        or output.get("benchmark_metrics") != "forbidden"
        or output.get("raw_training_targets_retained") is not False
        or output.get("champion_eligible") is not False
    ):
        raise SequenceTrainingError("campaign contract differs")


def configure_determinism(seed: int) -> None:
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
        for line in handle:
            value = line.strip()
            if value.startswith(">"):
                if name is not None:
                    if name in records:
                        raise SequenceTrainingError("duplicate FASTA name")
                    records[name] = "".join(pieces)
                name, pieces = value[1:], []
            elif name is None:
                raise SequenceTrainingError("FASTA sequence precedes name")
            else:
                pieces.append(value.upper())
    if name is not None:
        if name in records:
            raise SequenceTrainingError("duplicate FASTA name")
        records[name] = "".join(pieces)
    return records


def encode_sequences(sequences: Sequence[str]) -> torch.Tensor:
    values = torch.zeros((len(sequences), 4, 230), dtype=torch.float32)
    for row_index, sequence in enumerate(sequences):
        if len(sequence) != 230 or set(sequence) - set(CHANNEL_INDEX):
            raise SequenceTrainingError("sequence geometry or alphabet differs")
        for position, base in enumerate(sequence):
            values[row_index, CHANNEL_INDEX[base], position] = 1.0
    return values


def target_arrays(
    entries: Sequence[Mapping[str, str]],
    targets: Sequence[Mapping[str, str]],
) -> tuple[torch.Tensor, torch.Tensor]:
    by_key: dict[tuple[str, str], tuple[float, float]] = {}
    for row in targets:
        key = (row["element_id"], row["context_id"])
        if key in by_key or row["context_id"] not in CONTEXTS:
            raise SequenceTrainingError("duplicate or unexpected target row")
        try:
            values = (
                float(row["reference_activity"]),
                float(row["alternative_activity"]),
            )
        except ValueError as error:
            raise SequenceTrainingError("target activity is not numeric") from error
        if not all(math.isfinite(value) for value in values):
            raise SequenceTrainingError("target activity is not finite")
        by_key[key] = values
    reference = torch.empty((len(entries), len(CONTEXTS)), dtype=torch.float32)
    alternative = torch.empty_like(reference)
    for row_index, entry in enumerate(entries):
        for context_index, context in enumerate(CONTEXTS):
            key = (entry["element_id"], context)
            if key not in by_key:
                raise SequenceTrainingError("training element lacks target")
            reference[row_index, context_index], alternative[
                row_index, context_index
            ] = by_key[key]
    expected_keys = {
        (entry["element_id"], context)
        for entry in entries
        for context in CONTEXTS
    }
    if set(by_key) != expected_keys:
        raise SequenceTrainingError("target file contains held or unused rows")
    return reference, alternative


def fit_standardization(
    reference: torch.Tensor, alternative: torch.Tensor, indices: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    values = torch.cat((reference[indices], alternative[indices]), dim=0)
    mean = values.mean(dim=0)
    sd = values.std(dim=0, unbiased=False)
    if not torch.isfinite(mean).all() or not torch.isfinite(sd).all() or (sd <= 0).any():
        raise SequenceTrainingError("training-only target standardization differs")
    return mean, sd


def composite_loss_rows(
    predicted_reference: torch.Tensor,
    predicted_alternative: torch.Tensor,
    target_reference: torch.Tensor,
    target_alternative: torch.Tensor,
) -> torch.Tensor:
    reference = F.smooth_l1_loss(
        predicted_reference, target_reference, beta=1.0, reduction="none"
    )
    alternative = F.smooth_l1_loss(
        predicted_alternative, target_alternative, beta=1.0, reduction="none"
    )
    delta = F.smooth_l1_loss(
        predicted_alternative - predicted_reference,
        target_alternative - target_reference,
        beta=1.0,
        reduction="none",
    )
    return 0.25 * reference + 0.25 * alternative + delta


def train_epochs(
    *,
    model: torch.nn.Module,
    reference_sequences: torch.Tensor,
    alternative_sequences: torch.Tensor,
    target_reference: torch.Tensor,
    target_alternative: torch.Tensor,
    training_indices: torch.Tensor,
    epochs: int,
    seed: int,
    phase_code: int,
    device: torch.device,
    validation_indices: torch.Tensor | None = None,
    validation_blocks: Sequence[str] | None = None,
    minimum_epochs: int = 10,
    patience: int = 20,
) -> tuple[list[dict[str, float]], int]:
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=0.0003, weight_decay=0.0001
    )
    history: list[dict[str, float]] = []
    best_loss = math.inf
    stale = 0
    stop_epoch = epochs
    for epoch in range(1, epochs + 1):
        model.train()
        generator = torch.Generator(device="cpu")
        generator.manual_seed(seed * 1_000_003 + phase_code * 10_007 + epoch)
        order = training_indices[
            torch.randperm(len(training_indices), generator=generator)
        ]
        for start in range(0, len(order), 64):
            index = order[start : start + 64]
            reverse_mask = torch.rand(
                len(index), generator=generator
            ) < 0.5
            batch_reference, batch_alternative = joint_orientation_augmentation(
                reference_sequences[index], alternative_sequences[index], reverse_mask
            )
            batch_reference = batch_reference.to(device, non_blocking=True)
            batch_alternative = batch_alternative.to(device, non_blocking=True)
            observed_reference = target_reference[index].to(device, non_blocking=True)
            observed_alternative = target_alternative[index].to(
                device, non_blocking=True
            )
            optimizer.zero_grad(set_to_none=True)
            predicted_reference = model(batch_reference)
            predicted_alternative = model(batch_alternative)
            loss = composite_loss_rows(
                predicted_reference,
                predicted_alternative,
                observed_reference,
                observed_alternative,
            ).mean()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
        if validation_indices is not None:
            if validation_blocks is None:
                raise SequenceTrainingError("validation block identities are missing")
            validation = validation_block_loss(
                model=model,
                reference_sequences=reference_sequences,
                alternative_sequences=alternative_sequences,
                target_reference=target_reference,
                target_alternative=target_alternative,
                indices=validation_indices,
                block_ids=validation_blocks,
                device=device,
            )
            history.append({"epoch": float(epoch), **validation})
            if epoch >= minimum_epochs:
                if validation["block_mean_loss"] < best_loss - 1e-12:
                    best_loss = validation["block_mean_loss"]
                    stale = 0
                else:
                    stale += 1
                if stale >= patience:
                    stop_epoch = epoch
                    break
    return history, stop_epoch


def validation_block_loss(
    *,
    model: torch.nn.Module,
    reference_sequences: torch.Tensor,
    alternative_sequences: torch.Tensor,
    target_reference: torch.Tensor,
    target_alternative: torch.Tensor,
    indices: torch.Tensor,
    block_ids: Sequence[str],
    device: torch.device,
) -> dict[str, float]:
    model.eval()
    row_losses: list[torch.Tensor] = []
    with torch.inference_mode():
        for start in range(0, len(indices), 64):
            index = indices[start : start + 64]
            activities = allele_pair_activity(
                model,
                reference_sequences[index].to(device, non_blocking=True),
                alternative_sequences[index].to(device, non_blocking=True),
                reverse_complement_average=True,
            )
            row_losses.append(
                composite_loss_rows(
                    activities["reference_orientation_mean"],
                    activities["alternative_orientation_mean"],
                    target_reference[index].to(device, non_blocking=True),
                    target_alternative[index].to(device, non_blocking=True),
                )
                .mean(dim=1)
                .cpu()
            )
    losses = torch.cat(row_losses).numpy()
    if len(losses) != len(block_ids):
        raise SequenceTrainingError("validation loss row census differs")
    by_block: dict[str, list[float]] = {}
    for block, loss in zip(block_ids, losses, strict=True):
        by_block.setdefault(block, []).append(float(loss))
    block_means = np.asarray(
        [np.mean(by_block[block]) for block in sorted(by_block)], dtype=np.float64
    )
    if len(block_means) < 2 or not np.isfinite(block_means).all():
        raise SequenceTrainingError("validation block-loss distribution differs")
    return {
        "block_mean_loss": float(block_means.mean()),
        "block_loss_standard_error": float(
            block_means.std(ddof=1) / math.sqrt(len(block_means))
        ),
        "validation_blocks": float(len(block_means)),
    }


def select_epoch_one_standard_error(
    history: Sequence[Mapping[str, float]], minimum_epochs: int
) -> dict[str, Any]:
    candidates = [row for row in history if int(row["epoch"]) >= minimum_epochs]
    if not candidates:
        raise SequenceTrainingError("no eligible inner-validation epoch")
    leader = min(candidates, key=lambda row: (row["block_mean_loss"], row["epoch"]))
    threshold = leader["block_mean_loss"] + leader["block_loss_standard_error"]
    selected = min(
        (row for row in candidates if row["block_mean_loss"] <= threshold),
        key=lambda row: row["epoch"],
    )
    return {
        "selected_epoch_count": int(selected["epoch"]),
        "leader_epoch": int(leader["epoch"]),
        "leader_block_mean_optimization_loss": leader["block_mean_loss"],
        "leader_block_loss_standard_error": leader[
            "block_loss_standard_error"
        ],
        "one_standard_error_threshold": threshold,
        "selection_rule": "earliest_epoch_within_one_block_standard_error",
    }


def checkpoint_state(model: torch.nn.Module) -> dict[str, torch.Tensor]:
    state = {
        key: value.detach().cpu().contiguous()
        for key, value in model.state_dict().items()
    }
    if not state or any(not isinstance(value, torch.Tensor) for value in state.values()):
        raise SequenceTrainingError("tensor-only checkpoint state differs")
    return state


def predict_test_rows(
    *,
    model: torch.nn.Module,
    entries: Sequence[Mapping[str, str]],
    indices: torch.Tensor,
    reference_sequences: torch.Tensor,
    alternative_sequences: torch.Tensor,
    mean: torch.Tensor,
    sd: torch.Tensor,
    model_id: str,
    implementation_id: str,
    seed: int,
    outer_fold: str,
    checkpoint_hash: str,
    device: torch.device,
) -> list[dict[str, Any]]:
    model.eval()
    rows: list[dict[str, Any]] = []
    with torch.inference_mode():
        for start in range(0, len(indices), 64):
            batch_indices = indices[start : start + 64]
            activity = allele_pair_activity(
                model,
                reference_sequences[batch_indices].to(device, non_blocking=True),
                alternative_sequences[batch_indices].to(device, non_blocking=True),
                reverse_complement_average=True,
            )
            reference = (
                activity["reference_orientation_mean"].cpu() * sd + mean
            )
            alternative = (
                activity["alternative_orientation_mean"].cpu() * sd + mean
            )
            for local_index, global_index in enumerate(batch_indices.tolist()):
                entry = entries[global_index]
                for context_index, context in enumerate(CONTEXTS):
                    reference_value = float(reference[local_index, context_index])
                    alternative_value = float(
                        alternative[local_index, context_index]
                    )
                    rows.append(
                        {
                            "element_id": entry["element_id"],
                            "source_locus_group_id": entry[
                                "source_locus_group_id"
                            ],
                            "long_range_block_id": entry["long_range_block_id"],
                            "outer_fold": outer_fold,
                            "model_id": model_id,
                            "implementation_id": implementation_id,
                            "seed": seed,
                            "context_id": context,
                            "predicted_reference_activity": format(
                                reference_value, ".9g"
                            ),
                            "predicted_alternative_activity": format(
                                alternative_value, ".9g"
                            ),
                            "predicted_alt_minus_ref_activity": format(
                                alternative_value - reference_value, ".9g"
                            ),
                            "prediction_units": (
                                "mean_replicate_log2_RNA_over_DNA_activity"
                            ),
                            "checkpoint_sha256": checkpoint_hash,
                            "benchmark_metrics_calculated": "false",
                        }
                    )
    return rows


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise SequenceTrainingError("outer-fold output exists")
    config = load_config(arguments.config)
    validate_campaign(config)
    task_tree = arguments.frozen_task.resolve(strict=True)
    if (
        file_sha256(task_tree / "ARTIFACTS.json")
        != config["frozen_task"]["artifacts_sha256"]
    ):
        raise SequenceTrainingError("frozen task authority differs")
    verify_frozen_tree(task_tree)
    task_path = task_tree / config["frozen_task"]["task_spec_member"]
    manifest_path = task_tree / config["frozen_task"]["fixture_manifest_member"]
    fasta_path = task_tree / config["frozen_task"]["sequence_fasta_member"]
    split_path = task_tree / config["frozen_task"]["split_plan_member"]
    for path, expected in (
        (task_path, config["frozen_task"]["task_spec_sha256"]),
        (manifest_path, config["frozen_task"]["fixture_manifest_sha256"]),
        (fasta_path, config["frozen_task"]["sequence_fasta_sha256"]),
        (split_path, config["frozen_task"]["split_plan_sha256"]),
    ):
        if file_sha256(path) != expected:
            raise SequenceTrainingError("frozen task member hash differs")
    task = load_json(task_path)
    if (
        task.get("metrics_calculated")
        or task.get("training_executed")
        or task.get("predictions_generated")
        or task.get("artifact_contract", {}).get("adapter_metrics_forbidden")
        is not True
    ):
        raise SequenceTrainingError("frozen TaskSpec firewall differs")
    target_receipt = load_json(arguments.training_targets / "receipt.json")
    if (
        target_receipt.get("status") != "pass_outer_training_targets_only"
        or target_receipt.get("outer_test_fold") != arguments.outer_fold
        or target_receipt.get("held_test_target_rows_emitted") != 0
        or target_receipt.get("held_test_outcomes_visible_to_training_adapter")
        is not False
        or target_receipt.get("benchmark_metrics_calculated")
    ):
        raise SequenceTrainingError("training-only target receipt differs")
    target_path = arguments.training_targets / "training_targets.tsv"
    if file_sha256(target_path) != target_receipt["target_view_sha256"]:
        raise SequenceTrainingError("training-only target view hash differs")
    manifest = read_tsv(manifest_path)
    split_rows = read_tsv(split_path)
    split = next(
        (row for row in split_rows if row["outer_test_fold"] == arguments.outer_fold),
        None,
    )
    if split is None or split["held_test_outcomes_visible_to_adapter"] != "false":
        raise SequenceTrainingError("outer split contract differs")
    records = read_fasta(fasta_path)
    if len(manifest) != 1033 or len(records) != 4132:
        raise SequenceTrainingError("fixture census differs")
    entries = sorted(manifest, key=lambda row: row["element_id"])
    reference_sequences = encode_sequences(
        [records[f"{row['element_id']}|REF"] for row in entries]
    )
    alternative_sequences = encode_sequences(
        [records[f"{row['element_id']}|ALT"] for row in entries]
    )
    final_training_entries = [
        row for row in entries if row["outer_fold"] != arguments.outer_fold
    ]
    final_reference, final_alternative = target_arrays(
        final_training_entries, read_tsv(target_path)
    )
    position = {row["element_id"]: index for index, row in enumerate(entries)}
    target_position = {
        row["element_id"]: index
        for index, row in enumerate(final_training_entries)
    }
    expanded_reference = torch.full((len(entries), 2), torch.nan)
    expanded_alternative = torch.full((len(entries), 2), torch.nan)
    for element, local_index in target_position.items():
        global_index = position[element]
        expanded_reference[global_index] = final_reference[local_index]
        expanded_alternative[global_index] = final_alternative[local_index]
    inner_valid_fold = split["inner_validation_fold"]
    inner_train_folds = set(split["inner_training_folds"].split(";"))
    inner_train_indices = torch.tensor(
        [
            index
            for index, row in enumerate(entries)
            if row["outer_fold"] in inner_train_folds
        ],
        dtype=torch.long,
    )
    inner_valid_indices = torch.tensor(
        [
            index
            for index, row in enumerate(entries)
            if row["outer_fold"] == inner_valid_fold
        ],
        dtype=torch.long,
    )
    final_train_indices = torch.tensor(
        [
            index
            for index, row in enumerate(entries)
            if row["outer_fold"] != arguments.outer_fold
        ],
        dtype=torch.long,
    )
    test_indices = torch.tensor(
        [
            index
            for index, row in enumerate(entries)
            if row["outer_fold"] == arguments.outer_fold
        ],
        dtype=torch.long,
    )
    if (
        torch.isnan(expanded_reference[final_train_indices]).any()
        or torch.isnan(expanded_alternative[final_train_indices]).any()
        or not torch.isnan(expanded_reference[test_indices]).all()
        or not torch.isnan(expanded_alternative[test_indices]).all()
        or set(inner_train_indices.tolist()).intersection(test_indices.tolist())
        or set(inner_valid_indices.tolist()).intersection(test_indices.tolist())
    ):
        raise SequenceTrainingError("held-test target isolation differs")
    device = torch.device(arguments.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise SequenceTrainingError("CUDA device requested but unavailable")
    arguments.output.mkdir(mode=0o750)
    checkpoint_root = arguments.output / "checkpoints"
    receipt_root = arguments.output / "training_receipts"
    checkpoint_root.mkdir(mode=0o750)
    receipt_root.mkdir(mode=0o750)
    prediction_rows: list[dict[str, Any]] = []
    fold_index = FOLDS.index(arguments.outer_fold)
    checkpoint_hashes: dict[str, list[str]] = {model: [] for model in MODEL_IDS}
    started = time.monotonic()
    for model_index, model_id in enumerate(MODEL_IDS):
        implementation_id = task["architectures"][model_id]["implementation_id"]
        model_checkpoint_dir = checkpoint_root / model_id
        model_receipt_dir = receipt_root / model_id
        model_checkpoint_dir.mkdir(mode=0o750)
        model_receipt_dir.mkdir(mode=0o750)
        for seed in SEEDS:
            fit_started = time.monotonic()
            configure_determinism(seed)
            if device.type == "cuda":
                torch.cuda.reset_peak_memory_stats(device)
            inner_mean, inner_sd = fit_standardization(
                expanded_reference, expanded_alternative, inner_train_indices
            )
            inner_reference = (expanded_reference - inner_mean) / inner_sd
            inner_alternative = (expanded_alternative - inner_mean) / inner_sd
            model = build_model(model_id, seed=seed, dropout=0.1).to(device)
            history, stop_epoch = train_epochs(
                model=model,
                reference_sequences=reference_sequences,
                alternative_sequences=alternative_sequences,
                target_reference=inner_reference,
                target_alternative=inner_alternative,
                training_indices=inner_train_indices,
                validation_indices=inner_valid_indices,
                validation_blocks=[
                    entries[index]["long_range_block_id"]
                    for index in inner_valid_indices.tolist()
                ],
                epochs=200,
                minimum_epochs=10,
                patience=20,
                seed=seed,
                phase_code=fold_index * 10 + model_index,
                device=device,
            )
            selection = select_epoch_one_standard_error(history, 10)
            del model
            if device.type == "cuda":
                torch.cuda.empty_cache()
            final_mean, final_sd = fit_standardization(
                expanded_reference, expanded_alternative, final_train_indices
            )
            final_reference_z = (expanded_reference - final_mean) / final_sd
            final_alternative_z = (expanded_alternative - final_mean) / final_sd
            configure_determinism(seed)
            model = build_model(model_id, seed=seed, dropout=0.1).to(device)
            train_epochs(
                model=model,
                reference_sequences=reference_sequences,
                alternative_sequences=alternative_sequences,
                target_reference=final_reference_z,
                target_alternative=final_alternative_z,
                training_indices=final_train_indices,
                epochs=selection["selected_epoch_count"],
                seed=seed,
                phase_code=50 + fold_index * 10 + model_index,
                device=device,
            )
            checkpoint = model_checkpoint_dir / f"seed-{seed}.safetensors"
            save_file(
                checkpoint_state(model),
                checkpoint,
                metadata={
                    "model_id": model_id,
                    "implementation_id": implementation_id,
                    "outer_fold": arguments.outer_fold,
                    "seed": str(seed),
                    "selected_epoch_count": str(selection["selected_epoch_count"]),
                    "benchmark_metrics_calculated": "false",
                },
            )
            checkpoint_hash = file_sha256(checkpoint)
            checkpoint_hashes[model_id].append(checkpoint_hash)
            prediction_rows.extend(
                predict_test_rows(
                    model=model,
                    entries=entries,
                    indices=test_indices,
                    reference_sequences=reference_sequences,
                    alternative_sequences=alternative_sequences,
                    mean=final_mean,
                    sd=final_sd,
                    model_id=model_id,
                    implementation_id=implementation_id,
                    seed=seed,
                    outer_fold=arguments.outer_fold,
                    checkpoint_hash=checkpoint_hash,
                    device=device,
                )
            )
            fit_receipt = {
                "schema_version": "masld-bench-sequence-control-fit-receipt-v1",
                "status": "pass_prediction_only_fit",
                "dataset_id": "gse281364",
                "model_id": model_id,
                "implementation_id": implementation_id,
                "outer_test_fold": arguments.outer_fold,
                "seed": seed,
                "seed_is_genuine_training_initialization": True,
                "inner_training_elements": len(inner_train_indices),
                "inner_validation_elements": len(inner_valid_indices),
                "inner_validation_blocks": len(
                    {entries[index]["long_range_block_id"] for index in inner_valid_indices.tolist()}
                ),
                "final_refit_elements": len(final_train_indices),
                "held_test_elements": len(test_indices),
                "held_test_outcomes_visible_to_adapter": False,
                "inner_target_mean": inner_mean.tolist(),
                "inner_target_sd": inner_sd.tolist(),
                "final_refit_target_mean": final_mean.tolist(),
                "final_refit_target_sd": final_sd.tolist(),
                "inner_stop_epoch": stop_epoch,
                **selection,
                "optimizer": "AdamW",
                "learning_rate": 0.0003,
                "weight_decay": 0.0001,
                "batch_size": 64,
                "checkpoint_format": "safetensors_tensor_only",
                "checkpoint_sha256": checkpoint_hash,
                "prediction_rows": len(test_indices) * len(CONTEXTS),
                "prediction_units": (
                    "unstandardized_mean_replicate_log2_RNA_over_DNA_activity"
                ),
                "reverse_complement_averaged": True,
                "benchmark_metrics_calculated": False,
                "evaluator_invoked": False,
                "model_ranked": False,
                "fit_elapsed_seconds": time.monotonic() - fit_started,
                "peak_cuda_memory_bytes": (
                    int(torch.cuda.max_memory_allocated(device))
                    if device.type == "cuda"
                    else 0
                ),
            }
            (model_receipt_dir / f"seed-{seed}.json").write_text(
                json.dumps(fit_receipt, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            print(
                json.dumps(
                    {
                        "status": "fit_complete",
                        "outer_fold": arguments.outer_fold,
                        "model_id": model_id,
                        "seed": seed,
                        "selected_epoch_count": selection[
                            "selected_epoch_count"
                        ],
                        "elapsed_seconds": fit_receipt["fit_elapsed_seconds"],
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
            del model
            if device.type == "cuda":
                torch.cuda.empty_cache()
        if len(set(checkpoint_hashes[model_id])) != len(SEEDS):
            raise SequenceTrainingError(
                "checkpoint hashes are not distinct across genuine seeds"
            )
    expected_prediction_rows = len(test_indices) * len(CONTEXTS) * len(SEEDS) * len(MODEL_IDS)
    if len(prediction_rows) != expected_prediction_rows:
        raise SequenceTrainingError("outer-fold prediction row census differs")
    prediction_rows.sort(
        key=lambda row: (
            row["model_id"],
            row["seed"],
            row["element_id"],
            row["context_id"],
        )
    )
    fields = tuple(prediction_rows[0])
    write_tsv(arguments.output / "predictions.tsv", fields, prediction_rows)
    redacted_target_receipt = {
        key: value
        for key, value in target_receipt.items()
        if key
        not in {
            "target_values",
            "reference_activity",
            "alternative_activity",
        }
    }
    (arguments.output / "target_materialization_receipt.json").write_text(
        json.dumps(redacted_target_receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    receipt = {
        "schema_version": "masld-bench-sequence-control-outer-fold-bundle-v1",
        "status": "pass_prediction_only_outer_fold_bundle",
        "dataset_id": "gse281364",
        "outer_test_fold": arguments.outer_fold,
        "models": list(MODEL_IDS),
        "seeds": list(SEEDS),
        "fit_count": len(MODEL_IDS) * len(SEEDS),
        "checkpoint_count": len(MODEL_IDS) * len(SEEDS),
        "training_receipt_count": len(MODEL_IDS) * len(SEEDS),
        "prediction_rows": len(prediction_rows),
        "held_test_elements": len(test_indices),
        "held_test_outcomes_visible_to_adapter": False,
        "raw_training_targets_retained": False,
        "benchmark_metrics_calculated": False,
        "evaluator_invoked": False,
        "model_ranking_performed": False,
        "checkpoint_hashes_distinct_within_model_outer_fold": True,
        "device": str(device),
        "cuda_device_name": (
            torch.cuda.get_device_name(device) if device.type == "cuda" else None
        ),
        "elapsed_seconds": time.monotonic() - started,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True), flush=True)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--frozen-task", type=Path, required=True)
    parser.add_argument("--training-targets", type=Path, required=True)
    parser.add_argument("--outer-fold", choices=FOLDS, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
