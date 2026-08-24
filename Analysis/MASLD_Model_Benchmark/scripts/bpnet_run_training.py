#!/usr/bin/env python3
"""Train the exact pinned, bias-free BPNet matched control."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random
import shutil
import sys
from typing import Any

import numpy as np

try:
    from .bpnet_contract import (
        BPNetContractError,
        DEFAULT_SEED,
        INPUT_LENGTH,
        OUTPUT_LENGTH,
        architecture_parameters,
        prepare_narrowpeak,
        read_chrom_sizes,
        read_fold,
        single_task,
        validate_model_inputs,
    )
except ImportError:
    from bpnet_contract import (
        BPNetContractError,
        DEFAULT_SEED,
        INPUT_LENGTH,
        OUTPUT_LENGTH,
        architecture_parameters,
        prepare_narrowpeak,
        read_chrom_sizes,
        read_fold,
        single_task,
        validate_model_inputs,
    )


SOURCE_REVISION = "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def build_signature(arguments: argparse.Namespace, fold: dict[str, list[str]]) -> str:
    payload = {
        "schema_version": "masld-bench-bpnet-training-signature-v1",
        "source_revision": SOURCE_REVISION,
        "seed": arguments.seed,
        "input_length": INPUT_LENGTH,
        "output_length": OUTPUT_LENGTH,
        "fold": fold,
        "epochs": arguments.epochs,
        "batch_size": arguments.batch_size,
        "threads": arguments.threads,
        "max_jitter": arguments.max_jitter,
        "negative_ratio": arguments.negative_ratio,
        "learning_rate": arguments.learning_rate,
        "min_learning_rate": arguments.min_learning_rate,
        "early_stop": arguments.early_stop,
        "early_stop_min_delta": arguments.early_stop_min_delta,
        "reduce_lr_patience": arguments.reduce_lr_patience,
        "reduce_lr_factor": arguments.reduce_lr_factor,
        "genome": str(arguments.genome.resolve(strict=True)),
        "chrom_sizes": str(arguments.chrom_sizes.resolve(strict=True)),
        "bigwig": str(arguments.bigwig.resolve(strict=True)),
        "regions": str(arguments.regions.resolve(strict=True)),
        "nonpeaks": str(arguments.nonpeaks.resolve(strict=True)),
        "fold_path": str(arguments.fold.resolve(strict=True)),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _training_peak_frame(path: Path):
    import pandas as pd

    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=("chrom", "start", "end", "name", "score", "strand", "signal", "p", "q", "summit"),
    )
    if frame.empty:
        raise BPNetContractError("training peak frame is empty")
    frame["start"] = frame["start"] + frame["summit"] - OUTPUT_LENGTH // 2
    frame["end"] = frame["start"] + OUTPUT_LENGTH
    return frame[["chrom", "start", "end"]]


def _verify_architecture(model: Any) -> dict[str, object]:
    import tensorflow as tf

    if model.input_shape != (None, INPUT_LENGTH, 4):
        raise BPNetContractError(f"BPNet input geometry differs: {model.input_shape!r}")
    if [tuple(output.shape) for output in model.outputs] != [
        (None, OUTPUT_LENGTH, 1),
        (None, 1),
    ]:
        raise BPNetContractError("BPNet output geometry differs")
    if len(model.inputs) != 1 or model.inputs[0].name.split(":")[0] != "sequence":
        raise BPNetContractError("BPNet unexpectedly accepts a bias/control input")
    convolutions = [
        layer for layer in model.layers if isinstance(layer, tf.keras.layers.Conv1D)
    ]
    filters = [layer.filters for layer in convolutions]
    dilations = [layer.dilation_rate[0] for layer in convolutions]
    if filters != [64] * 9 + [1]:
        raise BPNetContractError(f"BPNet convolution filters differ: {filters!r}")
    if dilations != [1, 2, 4, 8, 16, 32, 64, 128, 256, 1]:
        raise BPNetContractError(f"BPNet dilation census differs: {dilations!r}")
    return {
        "parameter_count": model.count_params(),
        "conv_filters": filters,
        "conv_dilations": dilations,
        "bias_or_control_input": False,
    }


def run(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists() and arguments.resume_from is None:
        raise BPNetContractError(f"refusing to overwrite {arguments.output}")
    if arguments.output.exists() and any(arguments.output.iterdir()):
        raise BPNetContractError("BPNet output directory must be new or empty")
    if arguments.epochs < 1 or arguments.batch_size < 1 or arguments.threads < 1:
        raise BPNetContractError("epochs, batch size, and threads must be positive")
    if not 0 <= arguments.max_jitter <= 128:
        raise BPNetContractError("maximum jitter must be between 0 and 128 bp")
    if arguments.seed < 0 or arguments.early_stop < 1:
        raise BPNetContractError("seed and early-stop values are invalid")

    model_inputs = (
        arguments.genome,
        arguments.chrom_sizes,
        arguments.bigwig,
        arguments.regions,
        arguments.nonpeaks,
        arguments.fold,
    )
    validate_model_inputs(model_inputs)
    fold = read_fold(arguments.fold)
    signature = build_signature(arguments, fold)

    arguments.output.mkdir(mode=0o750, parents=False, exist_ok=True)
    prepared = arguments.output / "prepared"
    prepared.mkdir(mode=0o750)
    admitted_contigs = set(fold["train"] + fold["valid"])
    train_contigs = set(fold["train"])
    chrom_sizes = read_chrom_sizes(arguments.chrom_sizes)
    jittered_flank = INPUT_LENGTH // 2 + arguments.max_jitter
    role_flank_by_contig = {
        **{contig: jittered_flank for contig in train_contigs},
        **{contig: INPUT_LENGTH // 2 for contig in fold["valid"]},
    }
    peak_count = prepare_narrowpeak(
        arguments.regions,
        prepared / "peaks.train_valid.narrowPeak",
        allowed_contigs=admitted_contigs,
        chrom_sizes=chrom_sizes,
        required_flank_by_contig=role_flank_by_contig,
    )
    training_peak_count = prepare_narrowpeak(
        arguments.regions,
        prepared / "peaks.train.narrowPeak",
        allowed_contigs=train_contigs,
        chrom_sizes=chrom_sizes,
        required_flank=jittered_flank,
    )
    nonpeak_count = prepare_narrowpeak(
        arguments.nonpeaks,
        prepared / "nonpeaks.train_valid.narrowPeak",
        allowed_contigs=admitted_contigs,
        chrom_sizes=chrom_sizes,
        required_flank_by_contig=role_flank_by_contig,
    )

    import tensorflow as tf

    from bpnet.generators.generators import MBPNetSequenceGenerator
    from bpnet.model import arch
    from bpnet.train import training as upstream_training
    from bpnet.utils.counts_loss_weight import get_recommended_counts_loss_weight

    random.seed(arguments.seed)
    np.random.seed(arguments.seed)
    tf.random.set_seed(arguments.seed)
    counts_loss_weight = float(
        get_recommended_counts_loss_weight(
            [str(arguments.bigwig.resolve(strict=True))],
            _training_peak_frame(prepared / "peaks.train.narrowPeak"),
            alpha=1.0,
            orig_multi_loss=False,
        )
    )
    parameters = architecture_parameters(counts_loss_weight)
    tasks = single_task(
        bigwig=arguments.bigwig,
        peaks=prepared / "peaks.train_valid.narrowPeak",
        nonpeaks=prepared / "nonpeaks.train_valid.narrowPeak",
        negative_ratio=arguments.negative_ratio,
    )
    atomic_json(prepared / "bpnet_architecture.json", parameters)
    atomic_json(prepared / "bpnet_tasks.json", tasks)
    atomic_json(prepared / "bpnet_fold.json", fold)

    batch_common = {
        "input_seq_len": INPUT_LENGTH,
        "output_len": OUTPUT_LENGTH,
        "max_jitter": arguments.max_jitter,
        "rev_comp_aug": True,
        "shuffle": True,
    }
    train_parameters = dict(batch_common, mode="train")
    valid_parameters = dict(
        batch_common, max_jitter=0, rev_comp_aug=False, shuffle=False, mode="val"
    )
    task_path = prepared / "bpnet_tasks.json"
    train_generator = MBPNetSequenceGenerator(
        str(task_path),
        train_parameters,
        str(arguments.genome),
        str(arguments.chrom_sizes),
        chroms=fold["train"],
        num_threads=arguments.threads,
        batch_size=arguments.batch_size,
        epochs=arguments.epochs,
        foreground_weight=1.0,
        background_weight=0.0,
    )
    valid_generator = MBPNetSequenceGenerator(
        str(task_path),
        valid_parameters,
        str(arguments.genome),
        str(arguments.chrom_sizes),
        chroms=fold["valid"],
        num_threads=arguments.threads,
        batch_size=arguments.batch_size,
        epochs=arguments.epochs,
        foreground_weight=1.0,
        background_weight=0.0,
    )
    if train_generator.len() < 1 or valid_generator.len() < 1:
        raise BPNetContractError("BPNet training or validation generator is empty")

    random.seed(arguments.seed)
    np.random.seed(arguments.seed)
    tf.random.set_seed(arguments.seed)
    model = arch.BPNet(tasks, parameters, orig_multi_loss=False, name_prefix="main")
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=arguments.learning_rate),
        loss=None,
    )
    architecture = _verify_architecture(model)

    history: list[dict[str, object]] = []
    start_epoch = 0
    best_epoch = 0
    best_loss = float("inf")
    validation_losses: list[float] = []
    lr_losses: list[float] = []
    checkpoint = tf.train.Checkpoint(model=model, optimizer=model.optimizer)
    checkpoint_directory = arguments.output / "training_checkpoints"
    manager = tf.train.CheckpointManager(checkpoint, checkpoint_directory, max_to_keep=2)
    best_weights = arguments.output / "bpnet.best.weights.h5"
    resume_validation: dict[str, object] | None = None
    restore_status = None
    early_stopped = False

    if arguments.resume_from is not None:
        resume_root = arguments.resume_from.resolve(strict=True)
        state = json.loads((resume_root / "training_state.json").read_text(encoding="utf-8"))
        if state["training_signature"] != signature:
            raise BPNetContractError("resume signature differs")
        restore_path = str(resume_root / state["checkpoint_relative_path"])
        restore_status = checkpoint.restore(restore_path)
        restore_status.assert_existing_objects_matched()
        start_epoch = int(state["completed_epochs"])
        best_epoch = int(state["best_epoch"])
        best_loss = float(state["best_validation_loss"])
        history = json.loads((resume_root / "history.json").read_text(encoding="utf-8"))
        validation_losses = [float(row["validation_loss"]) for row in history]
        lr_losses = list(validation_losses[-arguments.reduce_lr_patience - 1 :])
        shutil.copy2(resume_root / "bpnet.best.weights.h5", best_weights)
        expected_iterations = int(state["optimizer_iterations"])
        observed_iterations = int(model.optimizer.iterations.numpy())
        if observed_iterations != expected_iterations:
            raise BPNetContractError("restored optimizer iteration differs")
        probe_path = resume_root / state["resume_probe_relative_path"]
        if sha256_file(probe_path) != state["resume_probe_sha256"]:
            raise BPNetContractError("resume prediction probe digest differs")
        with np.load(probe_path, allow_pickle=False) as handle:
            expected_profile = handle["profile"]
            expected_count = handle["count"]
        probe_rng = np.random.default_rng(arguments.seed + 999)
        probe_sequence = np.eye(4, dtype=np.float32)[
            probe_rng.integers(0, 4, size=(2, INPUT_LENGTH))
        ]
        observed_profile, observed_count = model.predict(probe_sequence, verbose=0)
        restoration_max_abs = max(
            float(np.max(np.abs(expected_profile - observed_profile))),
            float(np.max(np.abs(expected_count - observed_count))),
        )
        if restoration_max_abs > 1e-6:
            raise BPNetContractError(
                f"restored checkpoint prediction differs: {restoration_max_abs}"
            )
        resume_validation = {
            "resume_from": str(resume_root),
            "restored_completed_epochs": start_epoch,
            "restored_optimizer_iterations": observed_iterations,
            "restoration_prediction_max_abs": restoration_max_abs,
            "restoration_prediction_atol": 1e-6,
        }

    stop_epoch = arguments.epochs
    if arguments.attempt_epoch_limit is not None:
        if arguments.attempt_epoch_limit < 1:
            raise BPNetContractError("attempt epoch limit must be positive")
        stop_epoch = min(arguments.epochs, start_epoch + arguments.attempt_epoch_limit)
    for epoch in range(start_epoch, stop_epoch):
        epoch_seed = arguments.seed + epoch
        random.seed(epoch_seed)
        np.random.seed(epoch_seed)
        tf.random.set_seed(epoch_seed)
        train_generator.curr_epoch = epoch
        valid_generator.curr_epoch = epoch
        train_result = model.fit(
            train_generator.gen(),
            epochs=1,
            steps_per_epoch=train_generator.len(),
            verbose=2,
        )
        if restore_status is not None:
            restore_status.assert_consumed()
            restore_status = None
        validation_result = model.evaluate(
            valid_generator.gen(),
            steps=valid_generator.len(),
            return_dict=True,
            verbose=2,
        )
        validation_loss = float(validation_result["loss"])
        if not np.isfinite(validation_loss):
            raise BPNetContractError("validation loss is non-finite")
        validation_losses.append(validation_loss)
        lr_losses.append(validation_loss)
        if validation_loss < best_loss:
            best_loss = validation_loss
            best_epoch = epoch + 1
            model.save_weights(best_weights)

        current_lr = float(model.optimizer.learning_rate.numpy())
        new_lr = upstream_training.reduce_lr_on_plateau(
            lr_losses,
            current_lr,
            factor=arguments.reduce_lr_factor,
            patience=arguments.reduce_lr_patience,
            min_lr=arguments.min_learning_rate,
        )
        if new_lr != current_lr:
            model.optimizer.learning_rate.assign(new_lr)
            lr_losses = [validation_loss]
        row = {
            "epoch": epoch + 1,
            "epoch_seed": epoch_seed,
            "learning_rate": current_lr,
            "next_learning_rate": float(model.optimizer.learning_rate.numpy()),
            "training_loss": float(train_result.history["loss"][0]),
            "validation_loss": validation_loss,
            "validation_profile_loss": float(
                validation_result["profile_predictions_loss"]
            ),
            "validation_logcounts_loss": float(
                validation_result["logcounts_predictions_loss"]
            ),
        }
        history.append(row)
        checkpoint_path = manager.save(checkpoint_number=epoch + 1)
        probe_rng = np.random.default_rng(arguments.seed + 999)
        probe_sequence = np.eye(4, dtype=np.float32)[
            probe_rng.integers(0, 4, size=(2, INPUT_LENGTH))
        ]
        probe_profile, probe_count = model.predict(probe_sequence, verbose=0)
        resume_probe_path = arguments.output / f"resume_probe_epoch_{epoch + 1}.npz"
        np.savez_compressed(
            resume_probe_path,
            profile=np.asarray(probe_profile, dtype=np.float32),
            count=np.asarray(probe_count, dtype=np.float32),
        )
        atomic_json(arguments.output / "history.json", history)
        atomic_json(
            arguments.output / "training_state.json",
            {
                "schema_version": "masld-bench-bpnet-training-state-v1",
                "training_signature": signature,
                "completed_epochs": epoch + 1,
                "best_epoch": best_epoch,
                "best_validation_loss": best_loss,
                "checkpoint_relative_path": Path(checkpoint_path).resolve().relative_to(
                    arguments.output.resolve()
                ).as_posix(),
                "optimizer_iterations": int(model.optimizer.iterations.numpy()),
                "resume_probe_relative_path": resume_probe_path.resolve().relative_to(
                    arguments.output.resolve()
                ).as_posix(),
                "resume_probe_sha256": sha256_file(resume_probe_path),
                "seed": arguments.seed,
            },
        )
        if upstream_training.early_stopping_check(
            validation_losses,
            patience=arguments.early_stop,
            min_delta=arguments.early_stop_min_delta,
        ):
            early_stopped = True
            break

    if not best_weights.is_file() or not history:
        raise BPNetContractError("BPNet training produced no best checkpoint")
    model.load_weights(best_weights)
    inference_model = tf.keras.Model(
        inputs=model.inputs, outputs=model.outputs, name="bpnet_sequence_only"
    )
    inference_path = arguments.output / "bpnet.inference.h5"
    inference_model.save(inference_path, include_optimizer=False, save_format="h5")
    final_weights = arguments.output / "bpnet.weights.h5"
    model.save_weights(final_weights)

    receipt = {
        "schema_version": "masld-bench-bpnet-training-receipt-v1",
        "status": "pass",
        "model_id": "bpnet",
        "source_revision": SOURCE_REVISION,
        "training_signature": signature,
        "seed": arguments.seed,
        "completed_epochs": len(history),
        "requested_epochs": arguments.epochs,
        "campaign_complete": early_stopped or len(history) >= arguments.epochs,
        "terminal_reason": (
            "early_stopping"
            if early_stopped
            else "requested_epochs_complete"
            if len(history) >= arguments.epochs
            else "attempt_epoch_limit"
        ),
        "best_epoch": best_epoch,
        "best_validation_loss": best_loss,
        "counts_loss_weight": counts_loss_weight,
        "train_steps": train_generator.len(),
        "validation_steps": valid_generator.len(),
        "peak_count_train_valid": peak_count,
        "peak_count_train": training_peak_count,
        "nonpeak_count_train_valid": nonpeak_count,
        "input_length": INPUT_LENGTH,
        "output_length": OUTPUT_LENGTH,
        "max_jitter": arguments.max_jitter,
        "reverse_complement_augmentation": True,
        "negative_ratio": arguments.negative_ratio,
        "boundary_policy": {
            "training_required_flank_bp": jittered_flank,
            "validation_required_flank_bp": INPUT_LENGTH // 2,
            "input_half_width_bp": INPUT_LENGTH // 2,
            "max_jitter_bp": arguments.max_jitter,
            "role_specific": True,
            "applied_to_peaks_and_nonpeaks": True,
        },
        "training_contigs": fold["train"],
        "validation_contigs": fold["valid"],
        "test_contigs_exposed": False,
        "observed_held_donor_atac_exposed": False,
        "evaluator_outcomes_exposed": False,
        "resume_validation": resume_validation,
        "architecture": architecture,
        "inference_model_sha256": sha256_file(inference_path),
        "weights_sha256": sha256_file(final_weights),
        "best_weights_sha256": sha256_file(best_weights),
    }
    atomic_json(arguments.output / "training_receipt.json", receipt)
    atomic_json(
        arguments.output / "arguments.json",
        {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(arguments).items()
        },
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--bigwig", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--nonpeaks", type=Path, required=True)
    parser.add_argument("--fold", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume-from", type=Path)
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--attempt-epoch-limit", type=int)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--learning-rate", type=float, default=0.001)
    parser.add_argument("--min-learning-rate", type=float, default=1e-5)
    parser.add_argument("--max-jitter", type=int, default=128)
    parser.add_argument("--negative-ratio", type=float, default=0.1)
    parser.add_argument("--early-stop", type=int, default=5)
    parser.add_argument("--early-stop-min-delta", type=float, default=0.001)
    parser.add_argument("--reduce-lr-patience", type=int, default=2)
    parser.add_argument("--reduce-lr-factor", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    arguments = parser.parse_args()
    try:
        result = run(arguments)
    except BPNetContractError as error:
        raise SystemExit(str(error)) from error
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
