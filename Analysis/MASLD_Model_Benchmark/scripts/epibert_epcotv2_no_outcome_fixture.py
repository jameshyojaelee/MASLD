#!/usr/bin/env python3
"""Build native-shape observed-ATAC inputs without labels or outcomes."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any


FORBIDDEN_KEY_PARTS = ("label", "outcome", "truth", "response", "observed_y")


class NoOutcomeFixtureError(ValueError):
    """Raised when a synthetic fixture could expose or encode an outcome."""


def _array_sha256(array: Any) -> str:
    import numpy as np

    value = np.ascontiguousarray(array)
    return sha256(value.tobytes(order="C")).hexdigest()


def _validate_keys(keys: list[str]) -> None:
    lowered = [key.lower() for key in keys]
    if any(part in key for key in lowered for part in FORBIDDEN_KEY_PARTS):
        raise NoOutcomeFixtureError("fixture contains a forbidden outcome-like key")


def _epibert_inputs() -> dict[str, Any]:
    import numpy as np

    sequence_codes = np.arange(524_288, dtype=np.int64) % 4
    sequence = np.eye(4, dtype=np.float32)[sequence_codes][None, :, :]
    atac = ((np.arange(131_072, dtype=np.float32) % 257.0) / 257.0)[
        None, :, None
    ]
    atac_mask = np.zeros((1, 131_072, 1), dtype=np.bool_)
    mask_start = (131_072 - 384) // 2
    atac_mask[:, mask_start : mask_start + 384, :] = True
    atac[atac_mask] = 0.0
    motif_context = np.linspace(-1.0, 1.0, 693, dtype=np.float32)[None, :]
    return {
        "sequence_one_hot": sequence,
        "observed_atac_masked_query": atac,
        "atac_mask": atac_mask,
        "global_motif_context": motif_context,
    }


def _epcotv2_inputs() -> dict[str, Any]:
    import numpy as np

    tensor = np.zeros((1, 600, 5, 1600), dtype=np.float32)
    positions = np.arange(600 * 1600, dtype=np.int64).reshape(600, 1600)
    for channel in range(4):
        tensor[0, :, channel, :] = (positions % 4 == channel).astype(np.float32)
    tensor[0, :, 4, :] = (positions % 251).astype(np.float32) / 251.0
    masked_variant_bin = np.zeros((1, 600), dtype=np.bool_)
    masked_variant_bin[:, 300] = True
    tensor[0, 300, 4, :] = 0.0
    return {
        "sequence_plus_observed_atac": tensor,
        "masked_variant_atac_bin": masked_variant_bin,
    }


def _runtime_probe(model_id: str, arrays: dict[str, Any], runtime: str) -> dict[str, Any]:
    import numpy as np

    if runtime == "numpy":
        version = np.__version__
        device = "cpu"
        checksum = sum(float(np.asarray(value, dtype=np.float64).sum()) for value in arrays.values())
    elif runtime == "tensorflow":
        import tensorflow as tf

        devices = tf.config.list_physical_devices("GPU")
        device = "/GPU:0" if devices else "/CPU:0"
        with tf.device(device):
            tensors = [tf.convert_to_tensor(value) for value in arrays.values()]
            checksum = sum(float(tf.reduce_sum(tf.cast(value, tf.float64)).numpy()) for value in tensors)
        version = tf.__version__
    elif runtime == "torch":
        import torch

        device = "cuda" if torch.cuda.is_available() else "cpu"
        tensors = [torch.from_numpy(value).to(device) for value in arrays.values()]
        checksum = sum(float(value.to(dtype=torch.float64).sum().item()) for value in tensors)
        version = torch.__version__
    else:
        raise NoOutcomeFixtureError("runtime is unsupported")
    if not math.isfinite(checksum):
        raise NoOutcomeFixtureError("runtime checksum is non-finite")
    expected_version = {"epibert": "2.13.0", "epcotv2": "2.2.1"}[model_id]
    return {
        "runtime": runtime,
        "version": version,
        "expected_native_version": expected_version,
        "exact_native_version": version == expected_version,
        "device": device,
        "finite_tensor_checksum": format(checksum, ".17g"),
    }


def build(model_id: str, output: Path, runtime: str) -> dict[str, Any]:
    import numpy as np

    if output.exists() or model_id not in {"epibert", "epcotv2"}:
        raise NoOutcomeFixtureError("model or output contract is invalid")
    arrays = _epibert_inputs() if model_id == "epibert" else _epcotv2_inputs()
    _validate_keys(list(arrays))
    runtime_receipt = _runtime_probe(model_id, arrays, runtime)
    output.mkdir(parents=True, mode=0o750)
    input_path = output / "inputs.npz"
    np.savez_compressed(input_path, **arrays)
    with np.load(input_path, allow_pickle=False) as loaded:
        _validate_keys(list(loaded.files))
        if set(loaded.files) != set(arrays):
            raise NoOutcomeFixtureError("frozen input key set differs")
        array_contract = {
            key: {
                "shape": list(loaded[key].shape),
                "dtype": str(loaded[key].dtype),
                "sha256": _array_sha256(loaded[key]),
            }
            for key in sorted(loaded.files)
        }
    if model_id == "epibert":
        query_semantics = {
            "model_category": "observed_atac_comparator",
            "observed_atac_required_at_query": True,
            "query_atac_topology": "same_context_as_prediction_query",
            "masked_accessibility_endpoint": True,
            "scored_atac_span_removed_from_input": True,
            "mask_span_bp": 1536,
            "input_atac_resolution_bp": 4,
            "global_motif_context_required": True,
            "rna_conditioning": False,
        }
    else:
        query_semantics = {
            "model_category": "observed_atac_comparator",
            "observed_atac_required_at_query": True,
            "query_atac_topology": "same_context_as_prediction_query",
            "atac_is_output_head": False,
            "variant_alleles_share_identical_atac_query": True,
            "variant_containing_atac_bin_masked_in_this_fixture": True,
            "published_eqtl_classifier_present": False,
            "rna_conditioning": False,
        }
    receipt = {
        "schema_version": "masld-bench-observed-atac-no-outcome-fixture-v1",
        "status": "pass",
        "model_id": model_id,
        "synthetic_inputs": True,
        "project_data_read": False,
        "outcomes_read": False,
        "labels_present": False,
        "model_checkpoint_loaded": False,
        "model_forward_executed": False,
        "runtime_input_tensor_probe": runtime_receipt,
        "array_contract": array_contract,
        "query_semantics": query_semantics,
        "rna_conditioned_atac_eligible": False,
        "terminal_disposition": (
            "input_fixture_passed_checkpoint_restore_and_numeric_parity_pending"
            if model_id == "epibert"
            else "input_fixture_passed_but_execution_blocked_terms"
        ),
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=("epibert", "epcotv2"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--runtime", choices=("numpy", "tensorflow", "torch"), required=True
    )
    arguments = parser.parse_args()
    result = build(arguments.model, arguments.output, arguments.runtime)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
