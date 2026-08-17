"""Audited same-seed GPU tolerance measurement for the paper CL setting."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .execution import require_execution_ownership, verify_execution_record
from .firewall import validate_program_firewall


POLICY_SCHEMA = "masld-cl-gpu-tolerance-policy-v41"


def _load_policy(config: dict[str, Any], value: str | Path, *, verify_hashes=True):
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "gpu_tolerance_policy_v41.json"
    )
    if path != expected:
        raise ContractError("GPU tolerance requires its source-controlled V41 policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("report_only_no_selection_threshold") is not True
        or set(policy.get("replicas", {})) != {"c", "e"}
        or policy.get("failed_hardware_attempts", {}).get("used_for_tolerance") is not False
    ):
        raise ContractError("GPU tolerance policy identity differs")
    sources = {}
    for replica, replica_spec in policy["replicas"].items():
        for name in ("spec", "source_lock", "execution_record", "nvidia_smi", "update_manifest"):
            source = (path.parent / replica_spec[name]["path"]).resolve()
            if verify_hashes and sha256_path(source) != replica_spec[name]["sha256"]:
                raise ContractError(f"GPU tolerance source changed: {replica}.{name}")
            sources[f"{replica}.{name}"] = source
    return path, policy, sources


def _arguments_without_output(arguments: list[str]) -> list[str]:
    result, index = [], 0
    while index < len(arguments):
        if arguments[index] == "--output":
            index += 2
            continue
        result.append(arguments[index])
        index += 1
    return result


def _numeric_tolerance(left: np.ndarray, right: np.ndarray) -> dict[str, Any]:
    left = np.asarray(left)
    right = np.asarray(right)
    if left.shape != right.shape or left.dtype != right.dtype:
        raise ContractError("GPU tolerance arrays have different shape or dtype")
    delta = np.abs(left.astype(np.float64) - right.astype(np.float64))
    return {
        "shape": list(left.shape),
        "dtype": str(left.dtype),
        "bitwise_identical": bool(np.array_equal(left, right)),
        "nonzero_elements": int(np.count_nonzero(delta)),
        "maximum_absolute_difference": float(delta.max(initial=0.0)),
        "mean_absolute_difference": float(delta.mean()) if delta.size else 0.0,
        "root_mean_square_difference": float(np.sqrt(np.mean(np.square(delta)))) if delta.size else 0.0,
        "absolute_difference_quantiles": {
            str(value): float(observed)
            for value, observed in zip(
                [0.5, 0.9, 0.99, 0.999, 0.9999],
                np.quantile(delta, [0.5, 0.9, 0.99, 0.999, 0.9999])
                if delta.size else np.zeros(5),
            )
        },
        "within_absolute_1e_6": bool(np.allclose(left, right, rtol=0, atol=1e-6)),
    }


def _collect_tensors(value, prefix=""):
    import torch

    output = {}
    if torch.is_tensor(value):
        output[prefix] = value.detach().cpu()
    elif isinstance(value, dict):
        for key, child in value.items():
            output.update(_collect_tensors(child, f"{prefix}/{key}"))
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            output.update(_collect_tensors(child, f"{prefix}/{index}"))
    return output


def measure_gpu_tolerance(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    import torch

    validate_program_firewall(config)
    policy_path, policy, sources = _load_policy(config, policy_value)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    manifests, records, specs, locks = {}, {}, {}, {}
    for replica in ("c", "e"):
        with sources[f"{replica}.spec"].open() as handle:
            specs[replica] = json.load(handle)
        with sources[f"{replica}.source_lock"].open() as handle:
            locks[replica] = json.load(handle)
        with sources[f"{replica}.update_manifest"].open() as handle:
            manifests[replica] = json.load(handle)
        records[replica] = verify_execution_record(
            sources[f"{replica}.execution_record"], pipeline_root,
            config["_config_sha256"], allow_historical_source=True,
        )
        manifest_path = sources[f"{replica}.update_manifest"]
        output_dir = manifest_path.parent
        owned = [
            manifest_path, output_dir / "embedding_manifest.json",
            output_dir / "embedding_latent.npy", output_dir / "embedding_cells.tsv.gz",
            output_dir / "training_latent.npy", output_dir / "model" / "model.pt",
        ]
        require_execution_ownership(records[replica], owned, role=f"GPU tolerance replica {replica}")
        text = sources[f"{replica}.nvidia_smi"].read_text()
        expected = policy["replicas"][replica]
        if (
            records[replica].get("status") != "complete"
            or records[replica].get("slurm_job_id") != expected["expected_job_id"]
            or expected["expected_gpu"] not in text
            or f"Driver Version                            : {expected['expected_driver']}" not in text
        ):
            raise ContractError(f"GPU tolerance execution or hardware differs: {replica}")

    if _arguments_without_output(specs["c"]["arguments"]) != _arguments_without_output(
        specs["e"]["arguments"]
    ):
        raise ContractError("GPU tolerance training arguments differ beyond output")
    if locks["c"]["source_identity"] != locks["e"]["source_identity"]:
        raise ContractError("GPU tolerance source identities differ")
    if locks["c"]["input_artifacts"] != locks["e"]["input_artifacts"]:
        raise ContractError("GPU tolerance input artifacts differ")

    required = policy["required_setting"]
    for replica, manifest in manifests.items():
        observed = {
            "method": manifest.get("method"),
            "model_kind": manifest.get("model_kind"),
            "seed": manifest.get("seed"),
            "ewc_lambda": manifest.get("ewc_lambda"),
            "replay_fraction": manifest.get("replay_fraction"),
            "query_datasets": manifest.get("query_datasets"),
            "stopping_epoch": manifest.get("stopping_epoch"),
            "cells": manifest.get("embedding", {}).get("n_cells"),
            "latent_dimensions": manifest.get("embedding", {}).get("n_latent"),
        }
        if observed != required:
            raise ContractError(f"GPU tolerance setting differs: {replica}: {observed}")

    left_dir = sources["c.update_manifest"].parent
    right_dir = sources["e.update_manifest"].parent
    embedding = _numeric_tolerance(
        np.load(left_dir / "embedding_latent.npy", mmap_mode="r"),
        np.load(right_dir / "embedding_latent.npy", mmap_mode="r"),
    )
    training = _numeric_tolerance(
        np.load(left_dir / "training_latent.npy", mmap_mode="r"),
        np.load(right_dir / "training_latent.npy", mmap_mode="r"),
    )
    left_model = torch.load(left_dir / "model" / "model.pt", map_location="cpu")
    right_model = torch.load(right_dir / "model" / "model.pt", map_location="cpu")
    left_tensors, right_tensors = _collect_tensors(left_model), _collect_tensors(right_model)
    if set(left_tensors) != set(right_tensors):
        raise ContractError("GPU tolerance checkpoint tensor names differ")
    unequal, maximum = [], 0.0
    for name in sorted(left_tensors):
        if left_tensors[name].shape != right_tensors[name].shape:
            raise ContractError(f"GPU tolerance checkpoint tensor shape differs: {name}")
        if not torch.equal(left_tensors[name], right_tensors[name]):
            unequal.append(name)
            maximum = max(maximum, float(
                (left_tensors[name].double() - right_tensors[name].double()).abs().max()
            ))
    checkpoints = {
        "tensor_names_identical": True,
        "tensor_count": len(left_tensors),
        "bitwise_identical_tensors": len(unequal) == 0,
        "unequal_tensor_count": len(unequal),
        "maximum_absolute_tensor_difference": maximum,
        "checkpoint_file_hashes_identical": sha256_path(left_dir / "model" / "model.pt")
        == sha256_path(right_dir / "model" / "model.pt"),
        "checkpoint_file_hash_difference_reason": (
            "serialization bytes differ despite identical named tensors"
            if not unequal else "one or more named tensors differ"
        ),
    }
    cells_identical = sha256_path(left_dir / "embedding_cells.tsv.gz") == sha256_path(
        right_dir / "embedding_cells.tsv.gz"
    )
    result = {
        "schema_version": "masld-cl-gpu-tolerance-v41",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "same_source_identity": True,
        "same_input_artifacts": True,
        "same_arguments_except_output": True,
        "same_cell_roster_and_order": cells_identical,
        "hardware": {
            replica: {
                "gpu": policy["replicas"][replica]["expected_gpu"],
                "driver": policy["replicas"][replica]["expected_driver"],
                "job_id": policy["replicas"][replica]["expected_job_id"],
            }
            for replica in ("c", "e")
        },
        "full_embedding": embedding,
        "training_embedding": training,
        "checkpoint": checkpoints,
        "observed_gpu_tolerance": embedding["maximum_absolute_difference"],
        "report_only_no_selection_threshold": True,
        "failed_blackwell_attempts_excluded": policy["failed_hardware_attempts"],
    }
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "gpu_tolerance.json", result)
    return result

