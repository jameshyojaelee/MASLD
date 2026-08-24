#!/usr/bin/env python3
"""Run one frozen sequence control through the native profile data path."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import random
from typing import Any, Mapping


MODEL_IDS = ("sequence_cnn_control", "sequence_transformer_control")


class SequenceControlTrainingError(ValueError):
    """Raised when a sequence-control training request violates its contract."""


def _load_config(path: Path, model_id: str) -> Mapping[str, Any]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if raw.get("schema_version") != "masld-bench-sequence-controls-v1":
        raise SequenceControlTrainingError("sequence-control config schema differs")
    shared = raw.get("shared_contract")
    architectures = raw.get("architectures")
    if not isinstance(shared, dict) or not isinstance(architectures, dict):
        raise SequenceControlTrainingError("sequence-control config fields differ")
    if model_id not in MODEL_IDS or model_id not in architectures:
        raise SequenceControlTrainingError("sequence-control model is not admitted")
    if (
        shared.get("input_length") != 2114
        or shared.get("profile_output_length") != 1000
        or shared.get("outer_split") != "donor0_genomic0"
        or shared.get("from_scratch") is not True
        or shared.get("held_donor_atac_available_to_model") is not False
        or shared.get("test_outcomes_used") is not False
    ):
        raise SequenceControlTrainingError("shared sequence-control contract differs")
    return architectures[model_id]


def _read_params(path: Path) -> dict[str, str]:
    parameters: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 2 or not fields[0] or fields[0] in parameters:
            raise SequenceControlTrainingError("base parameter TSV differs")
        parameters[fields[0]] = fields[1]
    required = {
        "counts_loss_weight",
        "filters",
        "n_dil_layers",
        "inputlen",
        "outputlen",
        "negative_sampling_ratio",
        "max_jitter",
        "chr_fold_path",
    }
    if not required.issubset(parameters):
        raise SequenceControlTrainingError("base parameter TSV is incomplete")
    if parameters["inputlen"] != "2114" or parameters["outputlen"] != "1000":
        raise SequenceControlTrainingError("base sequence geometry differs")
    return parameters


def _control_parameters(model_id: str, architecture: Mapping[str, Any]) -> dict[str, str]:
    common = {
        "control_model_id": model_id,
        "control_blocks": str(architecture["blocks"]),
        "control_dropout": str(architecture["dropout"]),
        "control_stem_kernel_size": str(architecture["stem_kernel_size"]),
    }
    if model_id == "sequence_cnn_control":
        common.update(
            {
                "control_channels": str(architecture["channels"]),
                "control_kernel_size": str(architecture["kernel_size"]),
            }
        )
    else:
        common.update(
            {
                "control_attention_dropout": str(architecture["attention_dropout"]),
                "control_ffn_width": str(architecture["ffn_width"]),
                "control_heads": str(architecture["heads"]),
                "control_profile_kernel_size": str(
                    architecture["profile_kernel_size"]
                ),
                "control_token_stride": str(architecture["token_stride"]),
                "control_width": str(architecture["width"]),
            }
        )
    return common


def _write_effective_params(
    output_prefix: Path, base_params: Path, model_id: str, architecture: Mapping[str, Any]
) -> Path:
    parameters = _read_params(base_params)
    overlap = set(parameters) & set(_control_parameters(model_id, architecture))
    if overlap:
        raise SequenceControlTrainingError("base TSV contains reserved control fields")
    parameters.update(_control_parameters(model_id, architecture))
    target = Path(f"{output_prefix}.effective.params.tsv")
    if target.exists():
        raise SequenceControlTrainingError("effective parameter TSV already exists")
    with target.open("x", encoding="utf-8") as handle:
        for key in sorted(parameters):
            handle.write(f"{key}\t{parameters[key]}\n")
    return target


def main() -> None:
    from chrombpnet.training.utils import argmanager

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--control-config", type=Path, required=True)
    parser.add_argument("--control-model-id", choices=MODEL_IDS, required=True)
    control, remaining = parser.parse_known_args()
    import sys

    original = sys.argv
    try:
        sys.argv = [original[0], *remaining]
        arguments = argmanager.fetch_train_args()
    finally:
        sys.argv = original

    if arguments.seed < 0 or arguments.epochs < 1 or arguments.early_stop < 1:
        raise SequenceControlTrainingError("training seed or epoch contract differs")
    for value in (
        control.control_config,
        Path(arguments.genome),
        Path(arguments.bigwig),
        Path(arguments.peaks),
        Path(arguments.nonpeaks),
        Path(arguments.chr_fold_path),
        Path(arguments.params),
        Path(arguments.architecture_from_file),
    ):
        if value.is_symlink() or not value.is_file():
            raise SequenceControlTrainingError("training input is not a regular file")
    output_prefix = Path(arguments.output_prefix)
    if output_prefix.parent.is_symlink() or not output_prefix.parent.is_dir():
        raise SequenceControlTrainingError("training output parent differs")
    if any(Path(f"{output_prefix}{suffix}").exists() for suffix in (".h5", ".log")):
        raise SequenceControlTrainingError("training output already exists")

    architecture = _load_config(control.control_config, control.control_model_id)
    effective = _write_effective_params(
        output_prefix,
        Path(arguments.params),
        control.control_model_id,
        architecture,
    )
    arguments.params = str(effective)
    arguments.control_config = str(control.control_config.resolve(strict=True))
    arguments.control_model_id = control.control_model_id
    os.environ["PYTHONHASHSEED"] = str(arguments.seed)

    import numpy as np
    import tensorflow as tf

    random.seed(arguments.seed)
    np.random.seed(arguments.seed)
    tf.random.set_seed(arguments.seed)
    from chrombpnet.training import train

    train.main(arguments)

    receipt = {
        "schema_version": "masld-bench-sequence-control-training-receipt-v1",
        "status": "pass",
        "model_id": control.control_model_id,
        "seed": arguments.seed,
        "input_length": 2114,
        "output_length": 1000,
        "from_scratch": True,
        "train_outcomes_used": True,
        "validation_outcomes_used_for_early_stopping": True,
        "test_outcomes_used": False,
        "held_donor_atac_exposed": False,
        "effective_params": effective.name,
    }
    Path(f"{output_prefix}.training_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
