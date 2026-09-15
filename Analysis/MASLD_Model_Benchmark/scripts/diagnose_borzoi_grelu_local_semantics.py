#!/usr/bin/env python3
"""Test whether four explicit semantics explain the bounded Borzoi mismatch."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping

from scripts.probe_borzoi_grelu_local_numeric_parity import (
    BorzoiNumericParityError,
    build_fixture,
    build_historical_model,
    compare_tensors,
    digest,
    validate_contract as validate_base_contract,
    verify_tree,
)
from scripts.reconcile_borzoi_grelu_to_local_port_mapping import map_key


class BorzoiSemanticsDiagnosticError(RuntimeError):
    """Raised when a diagnostic binding or fail-closed condition differs."""


def validate_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("schema_version") != (
        "masld-bench-borzoi-grelu-local-semantics-diagnostic-v1"
    ):
        raise BorzoiSemanticsDiagnosticError("diagnostic schema differs")
    changes = contract["diagnostic_changes"]
    comparison = contract["comparison"]
    execution = contract["execution_contract"]
    claims = contract["claim_boundary"]
    if (
        contract["model_id"] != "borzoi_grelu_converted_port"
        or changes["batch_norm_eps"] != 0.001
        or changes["layer_norm_eps"] != 0.001
        or changes["gelu_approximation"] != "tanh"
        or changes["attention_position_training_length"] != 4096
        or changes["changes_apply_only_to_historical_candidate"] is not True
        or changes["checkpoint_tensors_changed"] is not False
        or changes["state_key_mapping_changed"] is not False
        or changes["input_fixture_changed"] is not False
        or comparison["absolute_tolerance"] != 1e-6
        or comparison["relative_tolerance"] != 1e-5
        or comparison["completion_records_pass_or_fail"] is not True
        or comparison["success_is_explanatory_only"] is not True
        or comparison["failure_retains_unresolved_semantic_differences"] is not True
        or execution["torch_load_weights_only"] is not True
        or execution["torch_load_unrestricted"] is not False
        or execution["checkpoint_export_allowed"] is not False
        or execution["prediction_export_allowed"] is not False
        or execution["biological_data_allowed"] is not False
        or execution["sealed_sources_allowed"] is not False
        or execution["global_frozen_census_mutation_allowed"] is not False
        or execution["local_port_execution_gate_may_open"] is not False
        or claims["native_borzoi_substitute"] is not False
        or claims["checkpoint_bound_source_recovered"] is not False
        or claims["native_numeric_parity_established"] is not False
        or claims["open_champion_eligible"] is not False
        or claims["external_champion_claim_allowed"] is not False
        or claims["universal_claim_allowed"] is not False
    ):
        raise BorzoiSemanticsDiagnosticError("diagnostic contract opened or drifted")


def fixed_training_length_central_mask(value, out_channels: int, training_length: int):
    import torch

    sequence_length = value.shape[-2]
    features = out_channels // 2
    rate = torch.exp(
        torch.log(
            torch.tensor(
                [training_length], device=value.device, dtype=value.dtype
            )
            + 1
        )
        / features
    )
    positions = torch.arange(
        -sequence_length + 1,
        sequence_length,
        device=value.device,
        dtype=value.dtype,
    )
    widths = rate ** torch.arange(
        1, features + 1, device=value.device, dtype=value.dtype
    ) - 1
    embeddings = widths[None, ...] > positions.abs()[..., None]
    signed = torch.sign(positions)[..., None] * embeddings
    return torch.cat((embeddings, signed), dim=-1)


def align_historical_semantics(historical, changes: Mapping[str, Any]) -> dict[str, int]:
    import torch

    counts = {
        "batch_norm_eps_changed": 0,
        "layer_norm_eps_changed": 0,
        "gelu_approximation_changed": 0,
        "attention_position_functions_changed": 0,
    }
    for module in historical.modules():
        if isinstance(module, torch.nn.BatchNorm1d):
            module.eps = changes["batch_norm_eps"]
            counts["batch_norm_eps_changed"] += 1
        elif isinstance(module, torch.nn.LayerNorm):
            module.eps = changes["layer_norm_eps"]
            counts["layer_norm_eps_changed"] += 1
        if isinstance(module, torch.nn.GELU):
            module.approximate = changes["gelu_approximation"]
            counts["gelu_approximation_changed"] += 1
    for block in historical.embedding.transformer_tower.blocks:
        block.mha.positional_embed = lambda value, out_channels, length=changes[
            "attention_position_training_length"
        ]: fixed_training_length_central_mask(value, out_channels, length)
        counts["attention_position_functions_changed"] += 1
    expected = {
        "batch_norm_eps_changed": 11,
        "layer_norm_eps_changed": 16,
        "gelu_approximation_changed": 12,
        "attention_position_functions_changed": 8,
    }
    if counts != expected:
        raise BorzoiSemanticsDiagnosticError(
            f"diagnostic module cardinality differs: {counts}"
        )
    return counts


def compute_stage_pairs(historical, local, input_tensor):
    import torch

    with torch.inference_mode():
        hist_conv, hist_y0, hist_y1 = historical.embedding.conv_tower(input_tensor)
        hist_transformer = historical.embedding.transformer_tower(hist_conv)
        hist_unet = historical.embedding.unet_tower(
            hist_transformer, [hist_y0, hist_y1]
        )
        hist_embedding = historical.embedding.act(
            historical.embedding.pointwise_conv(hist_unet)
        )
        hist_raw = historical.head(hist_embedding)

        local_stem = local.conv_dna(input_tensor)
        local_y1 = local.res_tower(local_stem)
        local_y0 = local.unet1(local_y1)
        local_conv = local._max_pool(local_y0)
        local_transformer = local.transformer(local_conv.permute(0, 2, 1)).permute(
            0, 2, 1
        )
        local_unet = local.upsampling_unet1(local_transformer) + local.horizontal_conv1(
            local_y0
        )
        local_unet = local.separable1(local_unet)
        local_unet = local.upsampling_unet0(local_unet) + local.horizontal_conv0(
            local_y1
        )
        local_unet = local.separable0(local_unet)
        local_cropped = local.crop(local_unet.permute(0, 2, 1)).permute(0, 2, 1)
        local_embedding = local.final_joined_convs(local_cropped)
        local_raw = local.human_head(local_embedding.float())
    return {
        "conv_tower_output": (hist_conv, local_conv),
        "conv_skip_y0": (hist_y0, local_y0),
        "conv_skip_y1": (hist_y1, local_y1),
        "transformer_output": (hist_transformer, local_transformer),
        "unet_output": (hist_unet, local_unet),
        "embedding_output": (hist_embedding, local_embedding),
        "raw_human_head": (hist_raw, local_raw),
        "softplus_human_head": (
            torch.nn.functional.softplus(hist_raw),
            torch.nn.functional.softplus(local_raw),
        ),
    }


def run_diagnostic(contract_path: Path, project_root: Path) -> dict[str, Any]:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    validate_contract(contract)
    failed_root = verify_tree(project_root, contract["failed_numeric_fixture"])
    failed = json.loads(
        (failed_root / "probe/numeric_parity_fixture.json").read_text()
    )
    if (
        failed["comparison"]["all_stages_passed"] is not False
        or failed["terminal_disposition"]
        != contract["failed_numeric_fixture"]["required_disposition"]
    ):
        raise BorzoiSemanticsDiagnosticError("bound numeric failure differs")
    base_path = project_root / contract["base_contract"]
    if base_path.is_symlink() or digest(base_path) != contract["base_contract_sha256"]:
        raise BorzoiSemanticsDiagnosticError("base numeric contract differs")
    base = json.loads(base_path.read_text(encoding="utf-8"))
    validate_base_contract(base)
    mapping_root = verify_tree(project_root, base["semantic_mapping_artifact"])
    source_root = verify_tree(project_root, base["source_artifact"])

    os.environ["HF_HUB_OFFLINE"] = "1"
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    import torch
    from borzoi_pytorch import Borzoi

    fixture = base["fixture"]
    torch.set_num_threads(fixture["threads"])
    torch.use_deterministic_algorithms(True)
    historical = build_historical_model(mapping_root / "historical_sources")
    changes = align_historical_semantics(historical, contract["diagnostic_changes"])
    local = Borzoi.from_hparams(
        return_center_bins_only=True,
        bins_to_return=fixture["expected_output_bins"],
        enable_mouse_head=False,
        flashed=False,
    )
    state = torch.load(
        source_root / "sources" / fixture["checkpoint_filename"],
        map_location="cpu",
        weights_only=True,
        mmap=True,
    )
    historical.load_state_dict(state, strict=True, assign=True)
    remapped = {map_key(key)[0]: value for key, value in state.items()}
    local.load_state_dict(remapped, strict=True, assign=True)
    historical.eval()
    local.eval()
    input_tensor, _ = build_fixture(base)
    pairs = compute_stage_pairs(historical, local, input_tensor)
    atol = contract["comparison"]["absolute_tolerance"]
    rtol = contract["comparison"]["relative_tolerance"]
    comparisons = {
        name: compare_tensors(left, right, atol=atol, rtol=rtol)
        for name, (left, right) in pairs.items()
    }
    explained = all(row["all_elements_close"] for row in comparisons.values())
    terminal = (
        "four_explicit_semantic_differences_explain_bounded_mismatch_but_do_not_"
        "establish_authoritative_or_native_execution"
        if explained
        else "four_explicit_semantic_differences_do_not_fully_explain_bounded_"
        "mismatch_and_execution_semantics_remain_unresolved"
    )
    return {
        "schema_version": contract["schema_version"],
        "model_id": contract["model_id"],
        "failed_numeric_fixture_sha256": digest(failed_root / "ARTIFACTS.json"),
        "diagnostic_change_counts": changes,
        "comparison": {
            "absolute_tolerance": atol,
            "relative_tolerance": rtol,
            "stages": comparisons,
            "all_stages_passed": explained,
        },
        "weights_only_load": True,
        "unrestricted_torch_load_executed": False,
        "checkpoint_tensors_changed": False,
        "checkpoint_exported": False,
        "prediction_exported": False,
        "biological_data_accessed": False,
        "sealed_sources_accessed": False,
        "diagnostic_alignment_is_authoritative_execution": False,
        "native_numeric_parity_established": False,
        "local_port_execution_gate_opened": False,
        "open_champion_eligible": False,
        "terminal_disposition": terminal,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise BorzoiSemanticsDiagnosticError("refusing to overwrite diagnostic")
    receipt = run_diagnostic(arguments.contract, arguments.project_root)
    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "semantics_diagnostic.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
