#!/usr/bin/env python3
"""Train one project-local scBasset fold with genomic-valid early stopping."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
import random
import sys
from typing import Any


_NATIVE_PATH = (
    Path(__file__).parents[1]
    / "src"
    / "masld_bench"
    / "adapters"
    / "scbasset_native.py"
)
_NATIVE_SPEC = importlib.util.spec_from_file_location(
    "masld_scbasset_native_standalone", _NATIVE_PATH
)
if _NATIVE_SPEC is None or _NATIVE_SPEC.loader is None:
    raise RuntimeError(f"cannot load scBasset native adapter: {_NATIVE_PATH}")
native = importlib.util.module_from_spec(_NATIVE_SPEC)
sys.modules[_NATIVE_SPEC.name] = native
_NATIVE_SPEC.loader.exec_module(native)
BOTTLENECK_DIMENSION = native.BOTTLENECK_DIMENSION
SOURCE_REVISION = native.SOURCE_REVISION
hdf5_tensor_manifest = native.hdf5_tensor_manifest
make_model = native.make_model
one_hot_base_codes = native.one_hot_base_codes
sha256_file = native.sha256_file


class ScBassetTrainingError(RuntimeError):
    """Raised when local scBasset training violates its fold contract."""


def _cell_count(path: Path, expected_folds: set[int]) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "cell_id", "donor_id", "lineage", "outer_fold"
        ):
            raise ScBassetTrainingError("training-cell fields differ")
        rows = [dict(row) for row in reader]
    if not rows or len({row["cell_id"] for row in rows}) != len(rows):
        raise ScBassetTrainingError("training-cell roster is empty or duplicated")
    if {int(row["outer_fold"]) for row in rows} != expected_folds:
        raise ScBassetTrainingError("training-cell donor folds differ")
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input-artifacts-sha256", required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--split-id", default="donor0_genomic0")
    parser.add_argument("--donor-test-fold", type=int, default=0)
    parser.add_argument("--donor-valid-fold", type=int, default=1)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--early-stopping-patience", type=int, default=10)
    arguments = parser.parse_args()

    import h5py
    import numpy as np
    from scipy import sparse
    import tensorflow as tf

    inputs = arguments.inputs.resolve(strict=True)
    source = arguments.source.resolve(strict=True)
    if (
        arguments.donor_test_fold not in range(5)
        or arguments.donor_valid_fold not in range(5)
        or arguments.donor_test_fold == arguments.donor_valid_fold
        or arguments.split_id
        != f"donor{arguments.donor_test_fold}_genomic{arguments.donor_test_fold}"
    ):
        raise ScBassetTrainingError("split arguments differ")
    training_folds = set(range(5)).difference(
        {arguments.donor_test_fold, arguments.donor_valid_fold}
    )
    if arguments.output.exists():
        raise ScBassetTrainingError(f"output already exists: {arguments.output}")
    for value, label in (
        (arguments.input_artifacts_sha256, "input artifact"),
        (arguments.source_artifacts_sha256, "source artifact"),
    ):
        if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
            raise ScBassetTrainingError(f"{label} hash is invalid")
    if (
        arguments.batch_size < 1
        or arguments.max_epochs < 1
        or arguments.early_stopping_patience < 1
        or not math.isfinite(arguments.learning_rate)
        or arguments.learning_rate <= 0
    ):
        raise ScBassetTrainingError("training hyperparameters are invalid")
    summary = json.loads((inputs / "summary.json").read_text(encoding="utf-8"))
    if (
        summary.get("status") != "pass"
        or summary.get("split_id") != arguments.split_id
        or summary.get("donor_train_folds") != sorted(training_folds)
        or summary.get("donor_valid_fold") != arguments.donor_valid_fold
        or summary.get("donor_test_fold") != arguments.donor_test_fold
        or summary.get("held_donor_atac_read") is not False
        or summary.get("genomic_test_atac_parsed_only_to_exclude") is not True
        or summary.get("genomic_test_atac_used_for_matrix") is not False
        or summary.get("genomic_test_atac_exported") is not False
        or summary.get("model_training_input_contains_genomic_test_atac") is not False
        or summary.get("intermediate_lineage_matrices_removed") is not True
        or summary.get("missing_evidence_encoded_as_zero") is not False
    ):
        raise ScBassetTrainingError("input leakage summary differs")
    n_cells = _cell_count(inputs / "training_cells.tsv", training_folds)
    train_matrix = sparse.load_npz(inputs / "m_train.npz").tocsr()
    valid_matrix = sparse.load_npz(inputs / "m_valid.npz").tocsr()
    with h5py.File(inputs / "train_seqs.h5", "r") as handle:
        train_shape = tuple(handle["X"].shape)
    with h5py.File(inputs / "valid_seqs.h5", "r") as handle:
        valid_shape = tuple(handle["X"].shape)
    if (
        train_matrix.shape != (train_shape[0], n_cells)
        or valid_matrix.shape != (valid_shape[0], n_cells)
        or train_shape[1:] != (1344,)
        or valid_shape[1:] != (1344,)
        or np.any(train_matrix.data != 1)
        or np.any(valid_matrix.data != 1)
    ):
        raise ScBassetTrainingError("sequence/target axes or binary labels differ")

    class SparseSequence(tf.keras.utils.Sequence):
        def __init__(
            self,
            sequence_path: Path,
            labels: Any,
            *,
            batch_size: int,
            shuffle: bool,
            seed: int,
        ) -> None:
            self.sequence_path = sequence_path
            self.labels = labels
            self.batch_size = batch_size
            self.shuffle = shuffle
            self.order = np.arange(labels.shape[0], dtype=np.int64)
            self.generator = np.random.default_rng(seed)
            if shuffle:
                self.generator.shuffle(self.order)

        def __len__(self) -> int:
            return math.ceil(len(self.order) / self.batch_size)

        def __getitem__(self, index: int) -> tuple[Any, Any]:
            selected = self.order[
                index * self.batch_size : (index + 1) * self.batch_size
            ]
            with h5py.File(self.sequence_path, "r") as handle:
                codes = np.asarray(handle["X"][np.sort(selected)], dtype=np.int8)
            if not np.array_equal(selected, np.sort(selected)):
                inverse = np.argsort(np.argsort(selected))
                codes = codes[inverse]
            sequence = one_hot_base_codes(codes)
            targets = self.labels[selected].toarray().astype(np.float32, copy=False)
            return sequence, targets

        def on_epoch_end(self) -> None:
            if self.shuffle:
                self.generator.shuffle(self.order)

    random.seed(arguments.seed)
    np.random.seed(arguments.seed)
    tf.random.set_seed(arguments.seed)
    tf.keras.backend.clear_session()
    model = make_model(source, n_cells)
    model.compile(
        loss=tf.keras.losses.BinaryCrossentropy(),
        optimizer=tf.keras.optimizers.Adam(
            learning_rate=arguments.learning_rate,
            beta_1=0.95,
            beta_2=0.9995,
        ),
    )
    arguments.output.mkdir(parents=True, mode=0o750)
    weights = arguments.output / "best_model.weights.h5"
    history_path = arguments.output / "training_history.tsv"
    callbacks = [
        tf.keras.callbacks.TerminateOnNaN(),
        tf.keras.callbacks.ModelCheckpoint(
            str(weights),
            monitor="val_loss",
            mode="min",
            save_best_only=True,
            save_weights_only=True,
        ),
        tf.keras.callbacks.EarlyStopping(
            monitor="val_loss",
            mode="min",
            min_delta=1e-6,
            patience=arguments.early_stopping_patience,
            restore_best_weights=True,
        ),
        tf.keras.callbacks.CSVLogger(str(history_path), separator="\t"),
    ]
    train_sequence = SparseSequence(
        inputs / "train_seqs.h5",
        train_matrix,
        batch_size=arguments.batch_size,
        shuffle=True,
        seed=arguments.seed,
    )
    valid_sequence = SparseSequence(
        inputs / "valid_seqs.h5",
        valid_matrix,
        batch_size=arguments.batch_size,
        shuffle=False,
        seed=arguments.seed,
    )
    history = model.fit(
        train_sequence,
        validation_data=valid_sequence,
        epochs=arguments.max_epochs,
        callbacks=callbacks,
        workers=0,
        use_multiprocessing=False,
        verbose=2,
    )
    if not weights.is_file() or not history.history.get("val_loss"):
        raise ScBassetTrainingError("training did not produce a selected checkpoint")
    val_losses = [float(value) for value in history.history["val_loss"]]
    if any(not math.isfinite(value) for value in val_losses):
        raise ScBassetTrainingError("validation loss is non-finite")
    best_epoch = int(np.argmin(val_losses)) + 1
    model.load_weights(weights)
    members = hdf5_tensor_manifest(weights)
    (arguments.output / "best_model.h5_members.json").write_text(
        json.dumps(members, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema_version": "masld-bench-scbasset-fold-model-v1",
        "status": "pass",
        "model_id": "scbasset",
        "source_revision": SOURCE_REVISION,
        "split_id": arguments.split_id,
        "donor_test_fold": arguments.donor_test_fold,
        "donor_valid_fold": arguments.donor_valid_fold,
        "seed": arguments.seed,
        "n_training_cells": n_cells,
        "bottleneck_dimension": BOTTLENECK_DIMENSION,
        "parameter_count": int(model.count_params()),
        "batch_size": arguments.batch_size,
        "learning_rate": arguments.learning_rate,
        "max_epochs": arguments.max_epochs,
        "early_stopping_patience": arguments.early_stopping_patience,
        "epochs_completed": len(val_losses),
        "best_epoch": best_epoch,
        "best_genomic_valid_binary_crossentropy": min(val_losses),
        "selection_objective": "minimum_binary_crossentropy_on_genomic_valid_ccres_using_nested_training_donor_cell_outputs",
        "input_artifacts_sha256": arguments.input_artifacts_sha256,
        "source_artifacts_sha256": arguments.source_artifacts_sha256,
        "weights_sha256": sha256_file(weights),
        "hdf5_tensor_members": len(members),
        "training_initialization": "from_scratch",
        "public_tutorial_weights_used": False,
        "held_donor_atac_used": False,
        "genomic_test_atac_used": False,
        "genomic_test_atac_available_to_training": False,
        "held_cell_embedding_available": False,
    }
    (arguments.output / "model_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
