"""Lock the control-only v7 adapter weight before outcome evaluation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path
from .execution import source_identity
from .firewall import validate_program_firewall


def _eligible_summary(rows: list[dict[str, Any]], config: dict[str, Any]) -> dict[str, Any]:
    all_lineage = [row for row in rows if row["model_kind"] == "all_lineage"]
    lineages = [row for row in rows if row["model_kind"] in config["lineages"]]
    if len(all_lineage) != 1 or len(lineages) != len(config["lineages"]):
        raise ContractError("adapter setting lacks the exact six-model roster")
    preservation = all(
        row["gates"]["reference_macro_f1"]
        and row["gates"]["reference_neighborhood"]
        and row["gates"]["uniform_offset_invariant"]
        for row in rows
    )
    protocol = all(
        row["gates"]["unsorted_harmony_no_material_worsening"]
        and row["gates"]["unsorted_architecture_no_material_worsening"]
        for row in rows
    )
    all_lineage_alignment = (
        all_lineage[0]["gates"]["harmony_improvement"]
        and all_lineage[0]["gates"]["architecture_improvement"]
    )
    improved_lineages = sum(
        row["gates"]["harmony_improvement"]
        and row["gates"]["architecture_improvement"]
        for row in lineages
    )
    eligible = (
        preservation and protocol and all_lineage_alignment
        and improved_lineages >= config["gates"]["minimum_improved_lineages"]
    )
    return {
        "eligible": bool(eligible),
        "preservation": bool(preservation),
        "protocol_no_material_worsening": bool(protocol),
        "all_lineage_alignment": bool(all_lineage_alignment),
        "improved_lineages": int(improved_lineages),
        "mean_candidate_shift": float(np.mean([
            row["bootstrap"]["harmony"]["candidate_shift"] for row in rows
        ])),
    }


def select_control_adapter_weight(
    config: dict[str, Any], policy_value: str | Path,
    result_values: list[str | Path], output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    policy_path = Path(policy_value).resolve()
    if policy_path != pipeline_root / "reference" / "control_adapter_policy_v7.json":
        raise ContractError("adapter selection requires the source-controlled v7 policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != "masld-cl-control-adapter-policy-v7"
        or policy.get("config_sha256") != config["_config_sha256"]
    ):
        raise ContractError("adapter selection policy is invalid")
    grid = list(map(float, policy["global_reference_weight_grid"]))
    expected = len(grid) * (len(config["lineages"]) + 1)
    if len(result_values) != expected:
        raise ContractError("adapter selection requires the complete weight-by-model matrix")
    rows, sources = [], []
    for value in result_values:
        path = Path(value).resolve()
        with path.open() as handle:
            row = json.load(handle)
        if (
            row.get("schema_version") != "masld-cl-control-adapter-pilot-v7"
            or row.get("config_sha256") != config["_config_sha256"]
            or row.get("control_only") is not True
            or row.get("outcomes_unlocked") is not False
        ):
            raise ContractError("adapter selection received an invalid control result")
        manifest_path = Path(row["sources"]["adapter_manifest"])
        if sha256_path(manifest_path) != row["sources"]["adapter_manifest_sha256"]:
            raise ContractError("adapter selection source manifest changed")
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        if manifest.get("adapter_source_identity") != source_identity(pipeline_root):
            raise ContractError("adapter selection source identity changed")
        rows.append(row)
        sources.append({"path": str(path), "sha256": sha256_path(path)})
    keys = [(float(row["global_reference_weight"]), row["model_kind"]) for row in rows]
    expected_keys = {(weight, kind) for weight in grid for kind in ["all_lineage", *config["lineages"]]}
    if len(set(keys)) != len(keys) or set(keys) != expected_keys:
        raise ContractError("adapter selection matrix is duplicated or incomplete")
    summaries = {}
    for weight in grid:
        summaries[str(weight)] = _eligible_summary(
            [row for row in rows if float(row["global_reference_weight"]) == weight], config
        )
    eligible = [weight for weight in grid if summaries[str(weight)]["eligible"]]
    if not eligible:
        raise ContractError("no v7 adapter weight passes the control-only gates")
    selected = min(eligible, key=lambda weight: (summaries[str(weight)]["mean_candidate_shift"], weight))
    result = {
        "schema_version": "masld-cl-control-adapter-selection-v7",
        "config_sha256": config["_config_sha256"],
        "selection_frozen": True,
        "outcomes_unlocked": True,
        "selected_global_reference_weight": selected,
        "summaries": summaries,
        "policy_realpath": str(policy_path),
        "policy_sha256": sha256_path(policy_path),
        "control_results": sources,
        "source_identity": source_identity(pipeline_root),
    }
    result["lock_sha256"] = hashlib.sha256(canonical_json_bytes(result)).hexdigest()
    write_json_exclusive(output_value, result)
    return result
