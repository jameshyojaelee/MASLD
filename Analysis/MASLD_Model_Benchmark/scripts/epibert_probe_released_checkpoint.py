#!/usr/bin/env python3
"""Strict restore and deterministic forward probe for one released EpiBERT member."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib
import json
import math
from pathlib import Path
from typing import Any


FORBIDDEN_KEY_PARTS = ("label", "outcome", "truth", "response", "observed_y")


class EpiBERTCheckpointProbeError(ValueError):
    """Raised when an executable EpiBERT checkpoint contract differs."""


def constructor_contract(model_kind: str) -> dict[str, Any]:
    if model_kind not in {"pretrained_atac", "fine_tuned_rampage"}:
        raise EpiBERTCheckpointProbeError("model kind differs")
    common = {
        "kernel_transformation": "relu_kernel_transformation",
        "dropout_rate": 0.20,
        "pointwise_dropout_rate": 0.10,
        "input_length": 524_288,
        "output_length": 4_096,
        "num_heads": 8,
        "numerical_stabilizer": 1.0e-7,
        "max_seq_length": 4_096,
        "norm": True,
        "BN_momentum": 0.90,
        "normalize": True,
        "seed": 21 if model_kind == "pretrained_atac" else 19,
        "num_transformer_layers": 8,
        "final_point_scale": 6,
        "filter_list_seq": [512, 640, 640, 768, 896, 1024],
        "filter_list_atac": [32, 64],
        "use_rot_emb": True,
        "final_output_length": 4_092 if model_kind == "pretrained_atac" else 896,
    }
    if model_kind == "fine_tuned_rampage":
        common["predict_atac"] = True
    return common


def _validate_fixture(arrays: dict[str, Any]) -> None:
    import numpy as np

    if set(arrays) != {
        "sequence_one_hot",
        "observed_atac_masked_query",
        "atac_mask",
        "global_motif_context",
    }:
        raise EpiBERTCheckpointProbeError("fixture key roster differs")
    if any(
        part in key.lower() for key in arrays for part in FORBIDDEN_KEY_PARTS
    ):
        raise EpiBERTCheckpointProbeError("fixture contains an outcome-like key")
    sequence = np.asarray(arrays["sequence_one_hot"])
    atac = np.asarray(arrays["observed_atac_masked_query"])
    mask = np.asarray(arrays["atac_mask"])
    motif = np.asarray(arrays["global_motif_context"])
    if (
        sequence.shape != (1, 524_288, 4)
        or sequence.dtype != np.float32
        or atac.shape != (1, 131_072, 1)
        or atac.dtype != np.float32
        or mask.shape != atac.shape
        or mask.dtype != np.bool_
        or motif.shape != (1, 693)
        or motif.dtype != np.float32
        or int(mask.sum()) != 384
        or not np.all(atac[mask] == 0)
        or not np.all(sequence.sum(axis=-1) == 1)
    ):
        raise EpiBERTCheckpointProbeError("fixture shape, dtype, or mask differs")
    positions = np.flatnonzero(mask.reshape(-1))
    if positions.tolist() != list(range((131_072 - 384) // 2, (131_072 + 384) // 2)):
        raise EpiBERTCheckpointProbeError("fixture mask is not the central 1536-bp span")


def _array_sha256(array: Any) -> str:
    import numpy as np

    value = np.ascontiguousarray(array)
    return sha256(value.tobytes(order="C")).hexdigest()


def _motif_model_input(array: Any) -> Any:
    """Match the released deserializer's pre-batch ``(1, 693)`` tensor."""
    import numpy as np

    value = np.asarray(array)
    if value.shape != (1, 693) or value.dtype != np.float32:
        raise EpiBERTCheckpointProbeError("motif fixture contract differs")
    return value[:, np.newaxis, :]


def _output_contract(model_kind: str, output: Any) -> dict[str, Any]:
    import numpy as np

    if model_kind == "pretrained_atac":
        values = {"masked_atac": np.asarray(output, dtype=np.float32)}
        expected_shape = (1, 4_092, 1)
    else:
        if not isinstance(output, (tuple, list)) or len(output) != 2:
            raise EpiBERTCheckpointProbeError("fine-tuned output roster differs")
        values = {
            "masked_atac": np.asarray(output[0], dtype=np.float32),
            "rampage": np.asarray(output[1], dtype=np.float32),
        }
        expected_shape = (1, 896, 1)
    if any(
        value.shape != expected_shape
        or not np.all(np.isfinite(value))
        or np.any(value < 0)
        for value in values.values()
    ):
        raise EpiBERTCheckpointProbeError("model output shape or values differ")
    return {
        name: {
            "shape": list(value.shape),
            "dtype": str(value.dtype),
            "sha256": _array_sha256(value),
            "minimum": float(value.min()),
            "maximum": float(value.max()),
            "mean": float(value.mean()),
        }
        for name, value in values.items()
    }


def probe(
    model_kind: str,
    checkpoint_prefix: Path,
    fixture: Path,
    output: Path,
    seed: int,
) -> dict[str, Any]:
    import numpy as np
    import tensorflow as tf

    if output.exists() or seed < 0:
        raise EpiBERTCheckpointProbeError("output or seed contract differs")
    if not checkpoint_prefix.with_suffix(".index").is_file():
        raise EpiBERTCheckpointProbeError("checkpoint prefix differs")
    with np.load(fixture, allow_pickle=False) as loaded:
        arrays = {key: loaded[key] for key in loaded.files}
    _validate_fixture(arrays)
    devices = tf.config.list_physical_devices("GPU")
    if len(devices) != 1 or "L40S" not in tf.config.experimental.get_device_details(devices[0]).get("device_name", ""):
        raise EpiBERTCheckpointProbeError("exactly one L40S GPU is required")
    tf.keras.utils.set_random_seed(seed)
    tf.config.experimental.enable_op_determinism()
    tf.keras.mixed_precision.set_global_policy("mixed_bfloat16")
    module_name = (
        "src.models.epibert_atac_pretrain"
        if model_kind == "pretrained_atac"
        else "src.models.epibert_rampage_finetune"
    )
    model_module = importlib.import_module(module_name)
    model = model_module.epibert(**constructor_contract(model_kind))
    restore = tf.train.Checkpoint(model=model).restore(str(checkpoint_prefix))
    inputs = [
        tf.convert_to_tensor(arrays["sequence_one_hot"]),
        tf.convert_to_tensor(arrays["observed_atac_masked_query"]),
        tf.convert_to_tensor(_motif_model_input(arrays["global_motif_context"])),
    ]
    first = model(inputs, training=False)
    restore.assert_existing_objects_matched()
    restore.expect_partial()
    first_contract = _output_contract(model_kind, first)
    second = model(inputs, training=False)
    second_contract = _output_contract(model_kind, second)
    if first_contract != second_contract:
        raise EpiBERTCheckpointProbeError("repeated inference is not bit-identical")
    variable_elements = sum(math.prod(variable.shape) for variable in model.variables)
    trainable_elements = sum(
        math.prod(variable.shape) for variable in model.trainable_variables
    )
    if variable_elements <= 0 or trainable_elements <= 0:
        raise EpiBERTCheckpointProbeError("restored model variable inventory differs")
    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": "masld-bench-epibert-released-checkpoint-probe-v1",
        "status": "pass",
        "model_kind": model_kind,
        "checkpoint_prefix": str(checkpoint_prefix),
        "constructor_contract": constructor_contract(model_kind),
        "tensorflow": tf.__version__,
        "mixed_precision_policy": tf.keras.mixed_precision.global_policy().name,
        "gpu": tf.config.experimental.get_device_details(devices[0]).get(
            "device_name"
        ),
        "model_variable_count": len(model.variables),
        "model_variable_elements": variable_elements,
        "trainable_variable_count": len(model.trainable_variables),
        "trainable_variable_elements": trainable_elements,
        "restore_existing_objects_matched": True,
        "checkpoint_extra_training_state_expected": True,
        "published_command_num_heads": 4 if model_kind == "pretrained_atac" else 8,
        "released_checkpoint_inferred_num_heads": 8,
        "constructor_heads_bound_to_checkpoint_shape_audit": True,
        "outputs": first_contract,
        "repeat_inference_bit_identical": True,
        "backward_executed": False,
        "checkpoint_resume_executed": False,
        "synthetic_inputs": True,
        "project_data_read": False,
        "outcomes_read": False,
        "labels_read": False,
        "observed_atac_required_at_query": True,
        "motif_fixture_shape": [1, 693],
        "motif_model_input_shape": [1, 1, 693],
        "scored_atac_span_removed_from_input": True,
        "rna_conditioned_atac_eligible": False,
        "sealed_inference_eligible": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model-kind",
        choices=("pretrained_atac", "fine_tuned_rampage"),
        required=True,
    )
    parser.add_argument("--checkpoint-prefix", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    arguments = parser.parse_args()
    probe(
        arguments.model_kind,
        arguments.checkpoint_prefix,
        arguments.fixture,
        arguments.output,
        arguments.seed,
    )


if __name__ == "__main__":
    main()
