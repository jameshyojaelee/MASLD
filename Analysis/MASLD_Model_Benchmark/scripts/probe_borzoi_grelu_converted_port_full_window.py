#!/usr/bin/env python3
"""Run a synthetic full-window deterministic forward through the Borzoi converted port."""

from __future__ import annotations

import argparse
import csv
import gc
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

from scripts.reconcile_borzoi_grelu_to_local_port_mapping import map_key


class BorzoiFullWindowProbeError(RuntimeError):
    """Raised when a full-window probe binding or fail-closed condition differs."""


SCHEMA = "masld-bench-borzoi-grelu-converted-port-full-window-probe-v1"
TARGET_FIELDS = [
    "",
    "identifier",
    "file",
    "clip",
    "clip_soft",
    "scale",
    "sum_stat",
    "strand_pair",
    "description",
]
TARGET_INDICES = [22, 23, 868, 869, 878, 879, 1302, 1303, 1367, 1510, 1724, 2035]


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def tensor_digest(value: Any) -> str:
    array = value.detach().cpu().contiguous().numpy()
    return sha256(memoryview(array).cast("B")).hexdigest()


def verify_tree(project_root: Path, binding: Mapping[str, Any], label: str) -> Path:
    root = project_root / binding["path"]
    if (
        root.is_symlink()
        or not (root / "COMPLETE").is_file()
        or digest(root / "ARTIFACTS.json") != binding["artifacts_sha256"]
    ):
        raise BorzoiFullWindowProbeError(f"{label} artifact differs")
    return root


def validate_contract(contract: Mapping[str, Any]) -> None:
    fixture = contract["fixture"]
    execution = contract["execution_contract"]
    claims = contract["claim_boundary"]
    targets = contract["target_subset"]
    if (
        contract.get("schema_version") != SCHEMA
        or contract.get("model_id") != "borzoi_grelu_converted_port"
        or fixture["replicate"] != 0
        or fixture["checkpoint_filename"] != "human_state_dict_rep0.h5"
        or fixture["sequence_bp"] != 524288
        or fixture["sequence_pattern"] != "ACGT_repeated"
        or fixture["batch_size"] != 1
        or fixture["dtype"] != "float32"
        or fixture["device"] != "cuda"
        or fixture["gpu_name"] != "NVIDIA L40S"
        or fixture["return_center_bins_only"] is not True
        or fixture["output_bins"] != 6144
        or fixture["repeat_for_determinism"] != 2
        or fixture["seed"] != 20260825
        or fixture["flashed"] is not False
        or [row["index"] for row in targets] != TARGET_INDICES
        or execution["expected_source_tensor_count"] != 215
        or execution["expected_native_human_tracks"] != 7611
        or execution["semantic_key_mapping_required"] is not True
        or execution["state_dict_order_mapping_allowed"] is not False
        or execution["torch_load_weights_only"] is not True
        or execution["torch_load_unrestricted"] is not False
        or execution["strict_restoration_required"] is not True
        or execution["checkpoint_export_allowed"] is not False
        or execution["full_output_tensor_export_allowed"] is not False
        or execution["output_digest_and_scalar_summary_allowed"] is not True
        or execution["biological_data_allowed"] is not False
        or execution["development_outcomes_allowed"] is not False
        or execution["sealed_sources_allowed"] is not False
        or execution["training_allowed"] is not False
        or execution["benchmark_metrics_allowed"] is not False
        or execution["global_frozen_census_mutation_allowed"] is not False
        or execution["dispatcher_only_gpu_submission"] is not True
        or claims["converted_port_runtime_gate_only"] is not True
        or claims["native_borzoi_substitute"] is not False
        or claims["checkpoint_bound_source_recovered"] is not False
        or claims["native_numeric_parity_established"] is not False
        or claims["native_weight_authority_established"] is not False
        or claims["open_champion_eligible"] is not False
        or claims["external_champion_claim_allowed"] is not False
        or claims["universal_claim_allowed"] is not False
    ):
        raise BorzoiFullWindowProbeError("full-window contract opened or drifted")


def verify_runtime(contract: Mapping[str, Any]) -> tuple[Path, Path]:
    runtime = contract["runtime"]
    environment = Path(runtime["environment"])
    if environment.is_symlink() or not environment.is_dir():
        raise BorzoiFullWindowProbeError("runtime environment differs")
    for key in (
        "python",
        "port_model_source",
        "port_config_source",
        "port_transformer_source",
        "target_manifest",
    ):
        path = environment / runtime[key]
        if path.is_symlink() or digest(path) != runtime[f"{key}_sha256"]:
            raise BorzoiFullWindowProbeError(f"runtime binding differs: {key}")
    return environment, environment / runtime["target_manifest"]


def read_target_subset(path: Path, expected: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    selected: dict[int, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header != TARGET_FIELDS:
            raise BorzoiFullWindowProbeError("target manifest header differs")
        count = 0
        for count, fields in enumerate(reader, start=1):
            if len(fields) != len(TARGET_FIELDS) or int(fields[0]) != count - 1:
                raise BorzoiFullWindowProbeError("target manifest geometry differs")
            index = int(fields[0])
            if index in TARGET_INDICES:
                selected[index] = {
                    "index": index,
                    "identifier": fields[1],
                    "strand_pair": int(fields[7]),
                    "description": fields[8],
                }
    if count != 7611:
        raise BorzoiFullWindowProbeError("native human target count differs")
    rows = [selected[index] for index in TARGET_INDICES]
    if rows != list(expected):
        raise BorzoiFullWindowProbeError("prespecified target subset differs")
    return rows


def build_synthetic_sequence(length: int, torch: Any) -> Any:
    positions = torch.arange(length, dtype=torch.int64)
    fixture = torch.zeros((1, 4, length), dtype=torch.float32)
    fixture[0, positions.remainder(4), positions] = 1.0
    if (
        tuple(fixture.shape) != (1, 4, 524288)
        or not bool(torch.all(fixture.sum(dim=1) == 1))
        or [int(value) for value in fixture.sum(dim=(0, 2)).tolist()]
        != [131072, 131072, 131072, 131072]
    ):
        raise BorzoiFullWindowProbeError("synthetic fixture differs")
    return fixture


def run_probe(contract_path: Path, project_root: Path) -> dict[str, Any]:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    validate_contract(contract)
    source = verify_tree(project_root, contract["source_artifact"], "source")
    verify_tree(project_root, contract["semantic_mapping_artifact"], "mapping")
    verify_tree(project_root, contract["semantics_diagnostic_artifact"], "diagnostic")
    _, target_manifest = verify_runtime(contract)
    targets = read_target_subset(target_manifest, contract["target_subset"])
    fixture_contract = contract["fixture"]
    checkpoint = source / "sources" / fixture_contract["checkpoint_filename"]
    if checkpoint.is_symlink() or digest(checkpoint) != fixture_contract["checkpoint_sha256"]:
        raise BorzoiFullWindowProbeError("checkpoint binding differs")

    os.environ["HF_HUB_OFFLINE"] = "1"
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    import torch
    from borzoi_pytorch import Borzoi

    if not torch.cuda.is_available():
        raise BorzoiFullWindowProbeError("CUDA is unavailable")
    if torch.cuda.get_device_name(0) != fixture_contract["gpu_name"]:
        raise BorzoiFullWindowProbeError("GPU identity differs")
    torch.manual_seed(fixture_contract["seed"])
    torch.cuda.manual_seed_all(fixture_contract["seed"])
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.use_deterministic_algorithms(True)
    model = Borzoi.from_hparams(
        return_center_bins_only=True,
        bins_to_return=fixture_contract["output_bins"],
        enable_mouse_head=False,
        flashed=False,
    )
    state = torch.load(
        checkpoint,
        map_location="cpu",
        weights_only=True,
        mmap=True,
    )
    if len(state) != contract["execution_contract"]["expected_source_tensor_count"]:
        raise BorzoiFullWindowProbeError("source tensor count differs")
    remapped = {}
    for source_key, tensor in state.items():
        target_key, _ = map_key(source_key)
        if target_key in remapped:
            raise BorzoiFullWindowProbeError("semantic mapping is not bijective")
        remapped[target_key] = tensor
    if set(remapped) != set(model.state_dict()):
        raise BorzoiFullWindowProbeError("semantic mapping is not onto the local port")
    model.load_state_dict(remapped, strict=True, assign=True)

    selected_indices = [row["index"] for row in targets]
    original_head = model.human_head
    subset_head = torch.nn.Conv1d(1920, len(selected_indices), 1)
    with torch.no_grad():
        subset_head.weight.copy_(original_head.weight[selected_indices])
        subset_head.bias.copy_(original_head.bias[selected_indices])
    model.human_head = subset_head
    del original_head, state, remapped
    gc.collect()
    model.eval().cuda()
    fixture = build_synthetic_sequence(fixture_contract["sequence_bp"], torch).cuda()
    torch.cuda.reset_peak_memory_stats()
    outputs = []
    with torch.inference_mode():
        for _ in range(fixture_contract["repeat_for_determinism"]):
            outputs.append(model(fixture).detach().cpu())
            torch.cuda.synchronize()
    expected_shape = (1, len(selected_indices), fixture_contract["output_bins"])
    if any(tuple(value.shape) != expected_shape for value in outputs):
        raise BorzoiFullWindowProbeError("full-window output geometry differs")
    if not bool(torch.isfinite(outputs[0]).all()) or not bool(torch.equal(*outputs)):
        raise BorzoiFullWindowProbeError("repeated full-window outputs differ")
    output = outputs[0]
    return {
        "schema_version": SCHEMA,
        "model_id": contract["model_id"],
        "status": "pass_synthetic_full_window_deterministic_converted_port_forward",
        "source_artifacts_sha256": contract["source_artifact"]["artifacts_sha256"],
        "mapping_artifacts_sha256": contract["semantic_mapping_artifact"]["artifacts_sha256"],
        "diagnostic_artifacts_sha256": contract["semantics_diagnostic_artifact"]["artifacts_sha256"],
        "checkpoint_sha256": fixture_contract["checkpoint_sha256"],
        "torch_version": torch.__version__,
        "cuda_version": torch.version.cuda,
        "gpu_name": torch.cuda.get_device_name(0),
        "sequence_bp": fixture_contract["sequence_bp"],
        "input_shape": list(fixture.shape),
        "output_shape": list(output.shape),
        "target_subset": targets,
        "target_manifest_sha256": contract["runtime"]["target_manifest_sha256"],
        "output_sha256": tensor_digest(output),
        "repeat_output_sha256": [tensor_digest(value) for value in outputs],
        "repeat_outputs_bit_identical": True,
        "output_min": float(output.min()),
        "output_max": float(output.max()),
        "output_mean": float(output.mean()),
        "peak_gpu_memory_bytes": int(torch.cuda.max_memory_allocated()),
        "weights_only_load": True,
        "unrestricted_torch_load_executed": False,
        "semantic_key_mapping_used": True,
        "strict_restoration_passed": True,
        "full_output_tensor_exported": False,
        "checkpoint_exported": False,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_sources_read": False,
        "training_performed": False,
        "benchmark_metrics_computed": False,
        "global_frozen_census_modified": False,
        "native_borzoi_substitute": False,
        "checkpoint_bound_source_recovered": False,
        "native_numeric_parity_established": False,
        "native_weight_authority_established": False,
        "open_champion_eligible": False,
        "external_champion_claim_allowed": False,
        "universal_claim_allowed": False,
        "terminal_disposition": (
            "converted_port_full_window_runtime_passed_but_authoritative_native_"
            "source_weight_parity_and_biological_task_gates_remain_closed"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise BorzoiFullWindowProbeError("refusing to overwrite probe output")
    receipt = run_probe(arguments.contract, arguments.project_root)
    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    printable = dict(receipt)
    printable["target_subset"] = len(receipt["target_subset"])
    print(json.dumps(printable, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
