"""Freeze the V8 control offset before matched-support outcome evaluation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path
from .execution import verify_source_identity_payload
from .firewall import validate_program_firewall


def _summary(rows: list[dict[str, Any]], config: dict[str, Any]) -> dict[str, Any]:
    all_lineage = [row for row in rows if row["model_kind"] == "all_lineage"]
    lineages = [row for row in rows if row["model_kind"] in config["lineages"]]
    if len(all_lineage) != 1 or len(lineages) != len(config["lineages"]):
        raise ContractError("V8 setting lacks the exact six-model roster")
    preservation = all(
        row["gates"]["reference_macro_f1"] and row["gates"]["reference_neighborhood"]
        and row["gates"]["uniform_offset_invariant"] for row in rows
    )
    protocol = all(
        row["gates"]["unsorted_harmony_no_material_worsening"]
        and row["gates"]["unsorted_architecture_no_material_worsening"] for row in rows
    )
    all_alignment = all_lineage[0]["gates"]["harmony_improvement"] and all_lineage[0]["gates"]["architecture_improvement"]
    improved_lineages = sum(
        row["gates"]["harmony_improvement"] and row["gates"]["architecture_improvement"]
        for row in lineages
    )
    eligible = preservation and protocol and all_alignment and improved_lineages >= config["gates"]["minimum_improved_lineages"]
    return {
        "eligible": bool(eligible), "preservation": bool(preservation),
        "protocol_no_material_worsening": bool(protocol),
        "all_lineage_alignment": bool(all_alignment),
        "improved_lineages": int(improved_lineages),
        "mean_candidate_shift": float(np.mean([row["bootstrap"]["harmony"]["candidate_shift"] for row in rows])),
    }


def select_orthogonal_bridge_weight(
    config: dict[str, Any], policy_value: str | Path,
    result_values: list[str | Path], output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    policy_path = Path(policy_value).resolve()
    if policy_path != pipeline_root / "reference" / "orthogonal_bridge_policy_v8.json":
        raise ContractError("V8 selection requires the source-controlled policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    grid = list(map(float, policy["control_offset_weight_grid"]))
    expected = len(grid) * (len(config["lineages"]) + 1)
    if len(result_values) != expected:
        raise ContractError("V8 selection requires the complete weight-by-model matrix")
    rows, sources = [], []
    frozen_source = None
    for value in result_values:
        path = Path(value).resolve()
        with path.open() as handle:
            row = json.load(handle)
        if (
            row.get("schema_version") != "masld-cl-orthogonal-bridge-pilot-v8"
            or row.get("config_sha256") != config["_config_sha256"]
            or row.get("control_only") is not True or row.get("outcomes_unlocked") is not False
        ):
            raise ContractError("V8 selection received an invalid control result")
        manifest_path = Path(row["sources"]["bridge_manifest"])
        if sha256_path(manifest_path) != row["sources"]["bridge_manifest_sha256"]:
            raise ContractError("V8 bridge manifest changed")
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        observed_source = verify_source_identity_payload(manifest.get("source_identity"))
        if frozen_source is None:
            frozen_source = observed_source
        elif observed_source != frozen_source:
            raise ContractError("V8 control candidates have mixed source identities")
        rows.append(row)
        sources.append({"path": str(path), "sha256": sha256_path(path)})
    keys = {(float(row["control_offset_weight"]), row["model_kind"]) for row in rows}
    expected_keys = {(weight, kind) for weight in grid for kind in ("all_lineage", *config["lineages"])}
    if len(keys) != len(rows) or keys != expected_keys:
        raise ContractError("V8 control grid is duplicated or incomplete")
    summaries = {
        str(weight): _summary([row for row in rows if float(row["control_offset_weight"]) == weight], config)
        for weight in grid
    }
    eligible = [weight for weight in grid if summaries[str(weight)]["eligible"]]
    if not eligible:
        raise ContractError("no V8 weight passes the control-only gates")
    selected = min(eligible, key=lambda weight: (summaries[str(weight)]["mean_candidate_shift"], weight))
    result = {
        "schema_version": "masld-cl-orthogonal-bridge-selection-v8",
        "config_sha256": config["_config_sha256"], "selection_frozen": True,
        "outcomes_unlocked": True, "selected_control_offset_weight": selected,
        "summaries": summaries, "policy_realpath": str(policy_path),
        "policy_sha256": sha256_path(policy_path), "control_results": sources,
        "source_identity": frozen_source,
    }
    result["lock_sha256"] = hashlib.sha256(canonical_json_bytes(result)).hexdigest()
    write_json_exclusive(output_value, result)
    return result
