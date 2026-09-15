#!/usr/bin/env python3
"""Predict fixed cCREs with scBasset and export sequence-only donor profiles."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
import random
import sys


ROOT = Path(__file__).parents[1]


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


native = _load_module(
    "masld_scbasset_native_predict",
    ROOT / "src" / "masld_bench" / "adapters" / "scbasset_native.py",
)
profile = _load_module(
    "masld_scbasset_profile_predict",
    ROOT / "src" / "masld_bench" / "adapters" / "scbasset_sequence_profile.py",
)


class ScBassetPredictionError(RuntimeError):
    """Raised when scBasset prediction does not meet its frozen requirements."""


def read_tsv(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise ScBassetPredictionError(f"{path.name} fields differ")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ScBassetPredictionError(f"{path.name} is empty")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--role", choices=("valid", "test"), required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--join-namespace", required=True)
    parser.add_argument("--input-artifacts-sha256", required=True)
    parser.add_argument("--model-artifacts-sha256", required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--split-id", default="donor0_genomic0")
    parser.add_argument("--donor-test-fold", type=int, default=0)
    parser.add_argument("--donor-valid-fold", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--profile-pseudocount", type=float, default=1e-8)
    arguments = parser.parse_args()

    import h5py
    import numpy as np
    import tensorflow as tf

    inputs = arguments.inputs.resolve(strict=True)
    model_root = arguments.model.resolve(strict=True)
    source = arguments.source.resolve(strict=True)
    if (
        arguments.donor_test_fold not in range(5)
        or arguments.donor_valid_fold not in range(5)
        or arguments.donor_test_fold == arguments.donor_valid_fold
        or arguments.split_id
        != f"donor{arguments.donor_test_fold}_genomic{arguments.donor_test_fold}"
    ):
        raise ScBassetPredictionError("split arguments differ")
    if arguments.output.exists() or arguments.batch_size < 1:
        raise ScBassetPredictionError("output or batch-size contract is invalid")
    if not math.isfinite(arguments.profile_pseudocount) or arguments.profile_pseudocount <= 0:
        raise ScBassetPredictionError("profile pseudocount is invalid")
    for value in (
        arguments.run_id,
        arguments.input_artifacts_sha256,
        arguments.model_artifacts_sha256,
        arguments.source_artifacts_sha256,
    ):
        if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
            raise ScBassetPredictionError("campaign or artifact hash is invalid")
    receipt = json.loads((model_root / "model_receipt.json").read_text(encoding="utf-8"))
    if (
        receipt.get("status") != "pass"
        or receipt.get("model_id") != "scbasset"
        or receipt.get("split_id") != arguments.split_id
        or receipt.get("donor_test_fold", 0) != arguments.donor_test_fold
        or receipt.get("donor_valid_fold", 1) != arguments.donor_valid_fold
        or receipt.get("public_tutorial_weights_used") is not False
        or receipt.get("held_donor_atac_used") is not False
        or receipt.get("genomic_test_atac_used") is not False
    ):
        raise ScBassetPredictionError("fitted model receipt differs")
    cells = read_tsv(inputs / "training_cells.tsv", profile.CELL_FIELDS)
    regions = read_tsv(inputs / f"{arguments.role}_regions.tsv", profile.REGION_FIELDS)
    donors = read_tsv(inputs / f"held_{arguments.role}_donors.tsv", profile.DONOR_FIELDS)
    if any(row["role"] != arguments.role for row in regions) or any(
        row["evaluation_role"] != arguments.role for row in donors
    ):
        raise ScBassetPredictionError("prediction role roster differs")
    sequence_path = inputs / f"{arguments.role}_seqs.h5"
    with h5py.File(sequence_path, "r") as handle:
        shape = tuple(handle["X"].shape)
    if shape != (len(regions), 1344) or len(cells) != int(receipt["n_training_cells"]):
        raise ScBassetPredictionError("prediction sequence or training-cell axis differs")

    seed = int(receipt["seed"])
    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)
    tf.keras.backend.clear_session()
    model = native.make_model(source, len(cells))
    weights = model_root / "best_model.weights.h5"
    if native.sha256_file(weights) != receipt["weights_sha256"]:
        raise ScBassetPredictionError("selected local HDF5 weights differ")
    model.load_weights(weights)
    arguments.output.mkdir(parents=True, mode=0o750)
    raw_path = arguments.output / f"{arguments.role}_training_cell_predictions.npy"
    raw = np.lib.format.open_memmap(
        raw_path,
        mode="w+",
        dtype=np.float32,
        shape=(len(regions), len(cells)),
    )
    with h5py.File(sequence_path, "r") as handle:
        codes = handle["X"]
        for start in range(0, len(regions), arguments.batch_size):
            end = min(start + arguments.batch_size, len(regions))
            sequence = native.one_hot_base_codes(codes[start:end])
            reverse_complement = sequence[:, ::-1, ::-1].copy()
            forward = np.asarray(model(sequence, training=False), dtype=np.float32)
            reverse = np.asarray(
                model(reverse_complement, training=False), dtype=np.float32
            )
            values = 0.5 * (forward + reverse)
            if (
                values.shape != (end - start, len(cells))
                or np.any(~np.isfinite(values))
                or np.any(values < 0)
                or np.any(values > 1)
            ):
                raise ScBassetPredictionError("scBasset RC-TTA output differs")
            raw[start:end] = values
    raw.flush()
    bundle = profile.export_prediction_bundle(
        sequence_predictions=np.load(raw_path, mmap_mode="r", allow_pickle=False),
        training_cells=cells,
        regions=regions,
        held_donors=donors,
        output=arguments.output / "bundle",
        run_id=arguments.run_id,
        namespace=arguments.join_namespace,
        donor_test_fold=arguments.donor_test_fold,
        donor_valid_fold=arguments.donor_valid_fold,
        model_artifact_sha256=arguments.model_artifacts_sha256,
        held_atac=None,
        held_cell_ids=None,
        pseudocount=arguments.profile_pseudocount,
    )
    result = {
        "schema_version": "masld-bench-scbasset-fold-predictions-v1",
        "status": "pass",
        "model_id": "scbasset",
        "split_id": arguments.split_id,
        "donor_test_fold": arguments.donor_test_fold,
        "donor_valid_fold": arguments.donor_valid_fold,
        "role": arguments.role,
        "regions": len(regions),
        "training_cells": len(cells),
        "held_donors": len(donors),
        "reverse_complement_tta": True,
        "training_cell_output_columns_only": True,
        "held_cell_embedding_available": False,
        "held_donor_atac_used": False,
        "rna_used": False,
        "public_tutorial_weights_used": False,
        "raw_predictions_sha256": native.sha256_file(raw_path),
        "raw_predictions_size_bytes": raw_path.stat().st_size,
        "prediction_bundle": bundle.relative_to(arguments.output).as_posix(),
        "input_artifacts_sha256": arguments.input_artifacts_sha256,
        "model_artifacts_sha256": arguments.model_artifacts_sha256,
        "source_artifacts_sha256": arguments.source_artifacts_sha256,
    }
    (arguments.output / "prediction_receipt.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
