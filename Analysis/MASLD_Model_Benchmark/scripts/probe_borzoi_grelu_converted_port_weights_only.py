#!/usr/bin/env python3
"""Probe weights-only loading and local-port schema compatibility for Borzoi."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from hashlib import sha256
import inspect
import json
import os
from pathlib import Path
from typing import Any, Mapping


class WeightsOnlyProbeError(RuntimeError):
    """Raised when the prospective safe-load requirement differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def canonical_hash(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def validate_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("schema_version") != (
        "masld-bench-borzoi-grelu-converted-port-weights-only-probe-v1"
    ):
        raise WeightsOnlyProbeError("probe schema differs")
    load = contract["load_contract"]
    claims = contract["claim_boundary"]
    if (
        contract["model_id"] != "borzoi_grelu_converted_port"
        or load["torch_load_weights_only"] is not True
        or load["torch_load_unrestricted"] is not False
        or load["map_location"] != "cpu"
        or load["model_forward_allowed"] is not False
        or load["checkpoint_export_allowed"] is not False
        or claims["native_borzoi_substitute"] is not False
        or claims["global_frozen_census_mutation_allowed"] is not False
        or claims["open_champion_eligible"] is not False
        or claims["external_champion_claim_allowed"] is not False
        or claims["universal_claim_allowed"] is not False
        or claims["sealed_sources_accessed"] is not False
        or claims["biological_outcomes_accessed"] is not False
    ):
        raise WeightsOnlyProbeError("fail-closed weights-only contract opened")
    if load["member_order"] != [
        "human_state_dict_rep0.h5",
        "human_state_dict_rep1.h5",
        "human_state_dict_rep2.h5",
        "human_state_dict_rep3.h5",
    ]:
        raise WeightsOnlyProbeError("checkpoint member order differs")


def validate_runtime(contract: Mapping[str, Any]) -> dict[str, Any]:
    runtime = contract["runtime"]
    environment = Path(runtime["environment"])
    if environment.is_symlink() or not environment.is_dir():
        raise WeightsOnlyProbeError("runtime environment differs")
    receipt: dict[str, Any] = {"environment": environment.as_posix(), "files": {}}
    for key in (
        "python",
        "conda_history",
        "port_metadata",
        "port_model_source",
        "port_config_source",
    ):
        path = environment / runtime[key]
        expected = runtime[f"{key}_sha256"]
        observed = digest(path)
        if path.is_symlink() or not path.is_file() or observed != expected:
            raise WeightsOnlyProbeError(f"runtime binding differs: {key}")
        receipt["files"][key] = {
            "path": path.as_posix(),
            "sha256": observed,
        }
    return receipt


def tensor_schema(state_dict: Mapping[str, Any], torch: Any) -> list[dict[str, Any]]:
    if not isinstance(state_dict, (dict, OrderedDict)):
        raise WeightsOnlyProbeError("weights-only payload is not a state dictionary")
    records = []
    for key, value in state_dict.items():
        if not isinstance(key, str) or not isinstance(value, torch.Tensor):
            raise WeightsOnlyProbeError("state dictionary contains a non-tensor value")
        records.append(
            {
                "key": key,
                "shape": list(value.shape),
                "dtype": str(value.dtype),
                "numel": int(value.numel()),
                "requires_grad": bool(value.requires_grad),
            }
        )
    return records


def probe(contract_path: Path, project_root: Path) -> dict[str, Any]:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    validate_contract(contract)
    runtime_receipt = validate_runtime(contract)
    source = project_root / contract["source_artifact"]["path"]
    manifest = source / "ARTIFACTS.json"
    if (
        source.is_symlink()
        or not (source / "COMPLETE").is_file()
        or digest(manifest) != contract["source_artifact"]["artifacts_sha256"]
    ):
        raise WeightsOnlyProbeError("static-audit source artifact differs")

    os.environ["HF_HUB_OFFLINE"] = "1"
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    import torch

    torch.set_num_threads(1)
    signature = inspect.signature(torch.load)
    supports_mmap = "mmap" in signature.parameters
    load_options = {
        "map_location": "cpu",
        "weights_only": True,
    }
    if supports_mmap and contract["load_contract"]["mmap_preferred"]:
        load_options["mmap"] = True

    schemas: list[list[dict[str, Any]]] = []
    members = []
    first_state_dict = None
    unsafe_scan = getattr(torch.serialization, "get_unsafe_globals_in_checkpoint", None)
    for filename in contract["load_contract"]["member_order"]:
        path = source / "sources" / filename
        if unsafe_scan is None:
            unsafe_globals: list[str] | str = "API_UNAVAILABLE"
        else:
            unsafe_globals = sorted(str(item) for item in unsafe_scan(path))
        state_dict = torch.load(path, **load_options)
        schema = tensor_schema(state_dict, torch)
        if len(schema) != contract["load_contract"]["expected_tensor_count"]:
            raise WeightsOnlyProbeError("weights-only tensor count differs")
        schema_sha256 = canonical_hash(schema)
        schemas.append(schema)
        members.append(
            {
                "filename": filename,
                "checkpoint_sha256": digest(path),
                "weights_only_load_passed": True,
                "unrestricted_torch_load_executed": False,
                "mmap_used": bool(load_options.get("mmap", False)),
                "unsafe_globals_static_api": unsafe_globals,
                "tensor_count": len(schema),
                "total_numel": sum(item["numel"] for item in schema),
                "state_schema_sha256": schema_sha256,
            }
        )
        if first_state_dict is None:
            first_state_dict = state_dict
        else:
            del state_dict
    schema_hashes = {member["state_schema_sha256"] for member in members}
    if len(schema_hashes) != 1:
        raise WeightsOnlyProbeError("replicate state schemas differ")

    from borzoi_pytorch import Borzoi

    model = Borzoi.from_hparams(
        return_center_bins_only=contract["load_contract"][
            "local_port_return_center_bins_only"
        ],
        enable_mouse_head=contract["load_contract"][
            "local_port_enable_mouse_head"
        ],
    )
    model_schema = tensor_schema(model.state_dict(), torch)
    model_schema_sha256 = canonical_hash(model_schema)
    checkpoint_schema = schemas[0]
    checkpoint_by_key = {item["key"]: item for item in checkpoint_schema}
    model_by_key = {item["key"]: item for item in model_schema}
    missing = sorted(set(model_by_key) - set(checkpoint_by_key))
    unexpected = sorted(set(checkpoint_by_key) - set(model_by_key))
    shape_mismatches = sorted(
        key
        for key in set(checkpoint_by_key) & set(model_by_key)
        if checkpoint_by_key[key]["shape"] != model_by_key[key]["shape"]
        or checkpoint_by_key[key]["dtype"] != model_by_key[key]["dtype"]
    )
    strict_restoration_passed = False
    strict_restoration_error = None
    if not missing and not unexpected and not shape_mismatches:
        try:
            model.load_state_dict(first_state_dict, strict=True, assign=True)
            strict_restoration_passed = True
        except (RuntimeError, TypeError) as error:
            strict_restoration_error = str(error)
    del model
    del first_state_dict

    return {
        "schema_version": (
            "masld-bench-borzoi-grelu-converted-port-weights-only-probe-v1"
        ),
        "model_id": contract["model_id"],
        "source_artifacts_sha256": contract["source_artifact"][
            "artifacts_sha256"
        ],
        "runtime": runtime_receipt,
        "python_executable": os.path.realpath(os.sys.executable),
        "torch_version": torch.__version__,
        "torch_load_signature": str(signature),
        "members": members,
        "replicate_schemas_identical": True,
        "checkpoint_state_schema_sha256": next(iter(schema_hashes)),
        "local_port_state_schema_sha256": model_schema_sha256,
        "local_port_tensor_count": len(model_schema),
        "strict_local_port_restoration_passed": strict_restoration_passed,
        "strict_local_port_missing_keys": missing,
        "strict_local_port_unexpected_keys": unexpected,
        "strict_local_port_shape_or_dtype_mismatches": shape_mismatches,
        "strict_local_port_restoration_error": strict_restoration_error,
        "unrestricted_torch_load_executed": False,
        "model_forward_executed": False,
        "checkpoint_exported": False,
        "native_borzoi_substitute": False,
        "open_champion_eligible": False,
        "external_champion_claim_allowed": False,
        "universal_claim_allowed": False,
        "terminal_disposition": (
            "weights_only_schema_probe_complete_local_port_compatibility_is_"
            + ("strict" if strict_restoration_passed else "not_strict")
            + "_native_parity_and_weight_authority_remain_required"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise WeightsOnlyProbeError("refusing to overwrite probe output")
    receipt = probe(args.contract, args.project_root)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "weights_only_probe.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
