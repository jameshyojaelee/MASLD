#!/usr/bin/env python3
"""Reconcile gReLU Borzoi state keys to the existing local PyTorch port."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import re
from typing import Any, Mapping


class BorzoiMappingError(RuntimeError):
    """Raised when the prospective semantic mapping is not exact."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def map_conv_suffix(suffix: str) -> str:
    if suffix.startswith("conv."):
        return "conv_layer." + suffix.removeprefix("conv.")
    if suffix.startswith("norm.layer."):
        return "norm." + suffix.removeprefix("norm.layer.")
    raise BorzoiMappingError(f"unmapped convolutional suffix: {suffix}")


def map_key(source_key: str) -> tuple[str, str]:
    match = re.fullmatch(r"embedding\.conv_tower\.blocks\.(\d+)\.(.+)", source_key)
    if match:
        block = int(match.group(1))
        suffix = match.group(2)
        prefixes = {
            0: "conv_dna",
            1: "res_tower.0",
            2: "res_tower.2",
            3: "res_tower.4",
            4: "res_tower.6",
            5: "res_tower.8",
            6: "unet1.1",
        }
        if block not in prefixes:
            raise BorzoiMappingError(f"unmapped convolutional block: {block}")
        if block == 0:
            if not suffix.startswith("conv."):
                raise BorzoiMappingError("stem contains an unexpected tensor")
            mapped = "conv_layer." + suffix.removeprefix("conv.")
        else:
            mapped = map_conv_suffix(suffix)
        return f"{prefixes[block]}.{mapped}", "convolutional_tower"

    match = re.fullmatch(
        r"embedding\.transformer_tower\.blocks\.(\d+)\.(.+)", source_key
    )
    if match:
        block = int(match.group(1))
        suffix = match.group(2)
        transformer_map = {
            "norm.layer.weight": "0.fn.0.weight",
            "norm.layer.bias": "0.fn.0.bias",
            "mha.rel_content_bias": "0.fn.1.rel_content_bias",
            "mha.rel_pos_bias": "0.fn.1.rel_pos_bias",
            "mha.to_k.weight": "0.fn.1.to_k.weight",
            "mha.to_out.weight": "0.fn.1.to_out.weight",
            "mha.to_out.bias": "0.fn.1.to_out.bias",
            "mha.to_q.weight": "0.fn.1.to_q.weight",
            "mha.to_pos_k.weight": "0.fn.1.to_rel_k.weight",
            "mha.to_v.weight": "0.fn.1.to_v.weight",
            "ffn.dense1.norm.layer.weight": "1.fn.0.weight",
            "ffn.dense1.norm.layer.bias": "1.fn.0.bias",
            "ffn.dense1.linear.weight": "1.fn.1.weight",
            "ffn.dense1.linear.bias": "1.fn.1.bias",
            "ffn.dense2.linear.weight": "1.fn.4.weight",
            "ffn.dense2.linear.bias": "1.fn.4.bias",
        }
        if suffix not in transformer_map:
            raise BorzoiMappingError(f"unmapped transformer tensor: {suffix}")
        family = (
            "transformer_feed_forward"
            if suffix.startswith("ffn.")
            else "transformer_attention"
        )
        return f"transformer.{block}.{transformer_map[suffix]}", family

    match = re.fullmatch(r"embedding\.unet_tower\.blocks\.(\d+)\.(.+)", source_key)
    if match:
        block = int(match.group(1))
        suffix = match.group(2)
        if block not in {0, 1}:
            raise BorzoiMappingError(f"unmapped U-net block: {block}")
        upsample = "upsampling_unet1.0" if block == 0 else "upsampling_unet0.0"
        horizontal = "horizontal_conv1" if block == 0 else "horizontal_conv0"
        separable = "separable1.conv_layer" if block == 0 else "separable0.conv_layer"
        if suffix.startswith("conv."):
            return f"{upsample}.{map_conv_suffix(suffix.removeprefix('conv.'))}", "unet_convolution"
        if suffix.startswith("channel_transform.conv.layer."):
            tail = suffix.removeprefix("channel_transform.conv.layer.")
            return f"{horizontal}.conv_layer.{tail}", "unet_channel_transform"
        if suffix.startswith("channel_transform.norm.layer."):
            tail = suffix.removeprefix("channel_transform.norm.layer.")
            return f"{horizontal}.norm.{tail}", "unet_channel_transform"
        if suffix.startswith("sconv.depthwise."):
            tail = suffix.removeprefix("sconv.depthwise.")
            return f"{separable}.0.{tail}", "unet_separable_convolution"
        if suffix.startswith("sconv.pointwise."):
            tail = suffix.removeprefix("sconv.pointwise.")
            return f"{separable}.1.{tail}", "unet_separable_convolution"
        raise BorzoiMappingError(f"unmapped U-net tensor: {suffix}")

    if source_key.startswith("embedding.pointwise_conv."):
        suffix = source_key.removeprefix("embedding.pointwise_conv.")
        return (
            f"final_joined_convs.0.{map_conv_suffix(suffix)}",
            "pointwise_projection",
        )
    if source_key.startswith("head.channel_transform.conv.layer."):
        suffix = source_key.removeprefix("head.channel_transform.conv.layer.")
        return f"human_head.{suffix}", "human_output_head"
    raise BorzoiMappingError(f"unmapped state key: {source_key}")


def validate_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("schema_version") != (
        "masld-bench-borzoi-grelu-to-local-port-mapping-probe-v1"
    ):
        raise BorzoiMappingError("mapping contract schema differs")
    mapping = contract["mapping_contract"]
    claims = contract["claim_boundary"]
    if (
        contract["model_id"] != "borzoi_grelu_converted_port"
        or contract["historical_source_candidate"][
            "checkpoint_producing_revision_proven"
        ]
        is not False
        or mapping["mapping_must_be_bijective"] is not True
        or mapping["shape_and_dtype_must_match"] is not True
        or mapping["state_dict_order_used_for_mapping"] is not False
        or mapping["torch_load_weights_only"] is not True
        or mapping["torch_load_unrestricted"] is not False
        or mapping["model_forward_allowed"] is not False
        or mapping["remapped_checkpoint_export_allowed"] is not False
        or claims["historical_source_is_checkpoint_bound"] is not False
        or claims["native_borzoi_substitute"] is not False
        or claims["global_frozen_census_mutation_allowed"] is not False
        or claims["open_champion_eligible"] is not False
        or claims["external_champion_claim_allowed"] is not False
        or claims["universal_claim_allowed"] is not False
        or claims["sealed_sources_accessed"] is not False
        or claims["biological_outcomes_accessed"] is not False
    ):
        raise BorzoiMappingError("fail-closed semantic-mapping contract opened")


def verify_tree(project_root: Path, binding: Mapping[str, Any]) -> Path:
    root = project_root / binding["path"]
    if (
        root.is_symlink()
        or not (root / "COMPLETE").is_file()
        or digest(root / "ARTIFACTS.json") != binding["artifacts_sha256"]
    ):
        raise BorzoiMappingError("bound artifact differs")
    return root


def reconcile(
    contract_path: Path, project_root: Path, historical_sources: Path
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    validate_contract(contract)
    source = verify_tree(project_root, contract["static_source_artifact"])
    verify_tree(project_root, contract["weights_only_probe"])
    history = contract["historical_source_candidate"]
    source_receipts = {}
    for name, expected in history["files"].items():
        path = historical_sources / name
        if (
            path.is_symlink()
            or not path.is_file()
            or path.stat().st_size != expected["size_bytes"]
            or digest(path) != expected["sha256"]
        ):
            raise BorzoiMappingError(f"historical source differs: {name}")
        source_receipts[name] = {
            "sha256": expected["sha256"],
            "size_bytes": expected["size_bytes"],
            "upstream_path": expected["path"],
        }

    runtime = contract["runtime"]
    environment = Path(runtime["environment"])
    for key in ("python", "port_model_source"):
        path = environment / runtime[key]
        if path.is_symlink() or digest(path) != runtime[f"{key}_sha256"]:
            raise BorzoiMappingError(f"runtime binding differs: {key}")
    os.environ["HF_HUB_OFFLINE"] = "1"
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    import torch
    from borzoi_pytorch import Borzoi

    torch.set_num_threads(1)
    model = Borzoi.from_hparams(return_center_bins_only=False, enable_mouse_head=False)
    target_state = model.state_dict()
    expected_count = contract["mapping_contract"]["expected_target_tensor_count"]
    if len(target_state) != expected_count:
        raise BorzoiMappingError("local target tensor count differs")

    member_receipts = []
    mapping_rows: list[dict[str, Any]] | None = None
    for replicate, filename in enumerate(
        [
            "human_state_dict_rep0.h5",
            "human_state_dict_rep1.h5",
            "human_state_dict_rep2.h5",
            "human_state_dict_rep3.h5",
        ]
    ):
        state = torch.load(
            source / "sources" / filename,
            map_location="cpu",
            weights_only=True,
            mmap=True,
        )
        if len(state) != contract["mapping_contract"]["expected_source_tensor_count"]:
            raise BorzoiMappingError("source tensor count differs")
        remapped = {}
        rows = []
        for source_key, tensor in state.items():
            target_key, family = map_key(source_key)
            if target_key in remapped:
                raise BorzoiMappingError(f"non-bijective target key: {target_key}")
            if target_key not in target_state:
                raise BorzoiMappingError(f"mapped target key absent: {target_key}")
            target = target_state[target_key]
            if tensor.shape != target.shape or tensor.dtype != target.dtype:
                raise BorzoiMappingError(f"mapped tensor schema differs: {source_key}")
            remapped[target_key] = tensor
            rows.append(
                {
                    "source_key": source_key,
                    "target_key": target_key,
                    "semantic_family": family,
                    "shape": list(tensor.shape),
                    "dtype": str(tensor.dtype),
                    "numel": int(tensor.numel()),
                }
            )
        if set(remapped) != set(target_state):
            missing = sorted(set(target_state) - set(remapped))
            raise BorzoiMappingError(f"mapping is not onto target state: {missing}")
        model.load_state_dict(remapped, strict=True, assign=True)
        mapping_sha256 = sha256(
            json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        if mapping_rows is None:
            mapping_rows = rows
        elif rows != mapping_rows:
            raise BorzoiMappingError("replicate semantic mappings differ")
        member_receipts.append(
            {
                "replicate": replicate,
                "filename": filename,
                "tensor_count": len(rows),
                "total_numel": sum(row["numel"] for row in rows),
                "mapping_sha256": mapping_sha256,
                "strict_load_passed": True,
                "weights_only_load": True,
                "unrestricted_torch_load_executed": False,
            }
        )
    assert mapping_rows is not None
    family_counts: dict[str, int] = {}
    for row in mapping_rows:
        family = row["semantic_family"]
        family_counts[family] = family_counts.get(family, 0) + 1
    if sorted(family_counts) != sorted(contract["mapping_contract"]["semantic_families"]):
        raise BorzoiMappingError("semantic mapping family roster differs")
    return (
        {
            "schema_version": "masld-bench-borzoi-grelu-to-local-port-mapping-probe-v1",
            "model_id": contract["model_id"],
            "historical_source_candidate": {
                "repository": history["repository"],
                "revision": history["revision"],
                "selection_rule": history["selection_rule"],
                "checkpoint_producing_revision_proven": False,
                "files": source_receipts,
            },
            "torch_version": torch.__version__,
            "members": member_receipts,
            "semantic_family_tensor_counts": dict(sorted(family_counts.items())),
            "mapping_is_bijective": True,
            "shape_and_dtype_match": True,
            "state_dict_order_used_for_mapping": False,
            "all_four_strict_loads_passed": True,
            "unrestricted_torch_load_executed": False,
            "model_forward_executed": False,
            "remapped_checkpoint_exported": False,
            "native_borzoi_substitute": False,
            "open_champion_eligible": False,
            "external_champion_claim_allowed": False,
            "universal_claim_allowed": False,
            "terminal_disposition": (
                "semantic_mapping_and_strict_weights_only_restoration_passed_"
                "but_historical_source_is_not_checkpoint_bound_and_numeric_native_"
                "parity_reference_track_and_weight_authority_gates_remain_closed"
            ),
        },
        mapping_rows,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--historical-sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise BorzoiMappingError("refusing to overwrite mapping output")
    receipt, rows = reconcile(
        args.contract, args.project_root, args.historical_sources
    )
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "mapping_probe.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with (args.output / "semantic_key_mapping.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        handle.write("source_key\ttarget_key\tsemantic_family\tshape\tdtype\tnumel\n")
        for row in rows:
            handle.write(
                "\t".join(
                    (
                        row["source_key"],
                        row["target_key"],
                        row["semantic_family"],
                        ",".join(str(value) for value in row["shape"]),
                        row["dtype"],
                        str(row["numel"]),
                    )
                )
                + "\n"
            )


if __name__ == "__main__":
    main()
