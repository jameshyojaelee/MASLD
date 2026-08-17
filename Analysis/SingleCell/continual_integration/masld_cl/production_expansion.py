"""Stage-blind six-model expansion of the confirmed V32 representation."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .control_adapter import _donor_balanced_centroid
from .embedding import load_embedding, matched_rows
from .firewall import validate_program_firewall
from .full_atlas_projection import _project_frozen_bridge


POLICY_SCHEMA = "masld-cl-production-expansion-policy-v35"


def _load_policy(config: dict[str, Any], value: str | Path):
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "production_expansion_policy_v35.json"
    )
    if path != expected:
        raise ContractError("production expansion requires the source-controlled V35 policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("adaptation", {}).get("minimum_donors") != 3
        or "query_control" not in policy.get("adaptation", {}).get("forbidden_inputs", [])
        or policy.get("firewall", {}).get("primary_reference_coordinates_remain_bitwise_frozen") is not True
        or policy.get("firewall", {}).get("descriptive_only_cells_are_unavailable_to_fit") is not True
    ):
        raise ContractError("production-expansion policy identity differs")
    model_specs = policy.get("models", {})
    if (
        set(model_specs) != {
            "all_lineage", "hepatocytes", "macrophages", "fibroblasts",
            "cholangiocytes", "t_cells",
        }
        or model_specs["all_lineage"].get("lineage_label") is not None
        or {spec.get("lineage_label") for name, spec in model_specs.items() if name != "all_lineage"}
        != set(config["lineages"])
    ):
        raise ContractError("production-expansion model roster differs")
    resolved = {}
    for name, source in policy["sources"].items():
        source_path = (path.parent / source["path"]).resolve()
        if sha256_path(source_path) != source["sha256"]:
            raise ContractError(f"production-expansion source changed: {name}")
        resolved[name] = source_path
    for model, spec in policy["models"].items():
        for key in ("v9_manifest", "v9_embedding"):
            source_path = (path.parent / spec[key]).resolve()
            if sha256_path(source_path) != spec[f"{key}_sha256"]:
                raise ContractError(f"production-expansion model source changed: {model}.{key}")
            resolved[f"{model}.{key}"] = source_path
    with resolved["v32_confirmation"].open() as handle:
        if json.load(handle).get("confirmation_pass") is not True:
            raise ContractError("production expansion requires a passing V32 confirmation")
    return path, policy, resolved


def _model_cells(full_cells, lineage_label: str | None):
    if lineage_label is None:
        return np.arange(len(full_cells), dtype=np.int64), full_cells.copy()
    keep = full_cells["audit_cell_type"].astype(str).to_numpy() == lineage_label
    positions = np.flatnonzero(keep)
    cells = full_cells.iloc[positions].reset_index(drop=True).copy()
    cells["row_index"] = np.arange(len(cells), dtype=np.int64)
    return positions, cells


def _select_stage_blind_offset(latent, cells, dataset: str, minimum: int):
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    eligible = cells["analysis_eligible"].to_numpy(dtype=bool)
    primary = cells["primary_query"].to_numpy(dtype=bool)
    datasets = cells["dataset"].astype(str).to_numpy()
    preparations = cells["preparation"].astype(str).to_numpy()
    query = eligible & ~reference & ~primary & (datasets == dataset)
    query_donors = cells.loc[query, "donor_id"].astype(str).nunique()
    if query_donors < minimum:
        return None
    compatible = []
    for preparation in sorted(set(preparations[query])):
        query_protocol = query & (preparations == preparation)
        reference_protocol = reference & (preparations == preparation)
        n_query = cells.loc[query_protocol, "donor_id"].astype(str).nunique()
        n_reference = cells.loc[reference_protocol, "donor_id"].astype(str).nunique()
        if n_query >= minimum and n_reference >= minimum:
            compatible.append((min(n_query, n_reference), preparation, query_protocol, reference_protocol))
    if compatible:
        _, preparation, fit_query, fit_reference = sorted(
            compatible, key=lambda value: (-value[0], value[1])
        )[0]
        fit_scope = f"compatible_protocol:{preparation}"
    else:
        n_reference = cells.loc[reference, "donor_id"].astype(str).nunique()
        if n_reference < minimum:
            return None
        fit_query, fit_reference = query, reference
        fit_scope = "all_preparations_fallback"
    query_centroid, n_query = _donor_balanced_centroid(latent, cells, fit_query)
    reference_centroid, n_reference = _donor_balanced_centroid(
        latent, cells, fit_reference
    )
    return {
        "dataset": dataset,
        "fit_scope": fit_scope,
        "n_query_donors": int(n_query),
        "n_reference_donors": int(n_reference),
        "offset": np.asarray(reference_centroid - query_centroid, dtype=np.float64),
        "apply_mask": (datasets == dataset) & ~reference & ~primary,
        "fit_mask": fit_query,
    }


def build_production_expansion(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy, sources = _load_policy(config, policy_value)
    _, raw, raw_cells = load_embedding(sources["full_raw_pca"])
    _, v33, full_cells = load_embedding(sources["full_v33"])
    if not np.array_equal(
        raw_cells["cell_id"].astype(str).to_numpy(),
        full_cells["cell_id"].astype(str).to_numpy(),
    ):
        raise ContractError("production expansion raw and routed full rosters differ")
    secondary_datasets = sorted(set(
        full_cells.loc[
            full_cells["analysis_eligible"]
            & ~full_cells["strict_reference"]
            & ~full_cells["primary_query"], "dataset"
        ].astype(str)
    ))
    output_root = Path(output_value).resolve()
    output_root.mkdir(parents=True, exist_ok=False)
    bundle_models = {}
    minimum = int(policy["adaptation"]["minimum_donors"])
    for model_kind, spec in policy["models"].items():
        raw_positions, cells = _model_cells(full_cells, spec["lineage_label"])
        observed = {
            "descriptive_cells": len(cells),
            "analyzed_cells": int(cells["analysis_eligible"].sum()),
            "projection_only_cells": int((~cells["analysis_eligible"]).sum()),
            "frozen_primary_reference_cells": int((cells["strict_reference"] | cells["primary_query"]).sum()),
        }
        expected = {key: spec[key] for key in observed}
        if observed != expected:
            raise ContractError(f"production-expansion {model_kind} census differs: {observed}")
        v9_manifest_path = sources[f"{model_kind}.v9_manifest"]
        with v9_manifest_path.open() as handle:
            v9_manifest = json.load(handle)
        v9_info, v9_latent, v9_cells = load_embedding(sources[f"{model_kind}.v9_embedding"])
        if (
            v9_manifest.get("schema_version") != "masld-cl-orthogonal-bridge-retarget-v9"
            or float(v9_manifest.get("global_reference_weight", -1)) != 0.75
            or v9_manifest.get("embedding") != v9_info
        ):
            raise ContractError(f"production-expansion V9 identity differs: {model_kind}")

        model_output = output_root / model_kind
        model_output.mkdir()
        latent_path = model_output / "embedding_latent.npy"
        latent = np.lib.format.open_memmap(
            latent_path, mode="w+", dtype=np.float32, shape=(len(cells), 30)
        )
        if model_kind == "all_lineage":
            latent[:] = np.asarray(v33, dtype=np.float32)
        else:
            with Path(v9_manifest["parent_manifest_realpath"]).open() as handle:
                bridge = json.load(handle)
            if sha256_path(v9_manifest["parent_manifest_realpath"]) != v9_manifest["parent_manifest_sha256"]:
                raise ContractError(f"production-expansion V8 parent changed: {model_kind}")
            _project_frozen_bridge(raw[raw_positions], bridge, latent)
            v9_left, full_right = matched_rows(v9_cells, cells)
            order = np.argsort(full_right)
            v9_left, full_right = v9_left[order], full_right[order]
            if len(full_right) != spec["frozen_primary_reference_cells"]:
                raise ContractError(f"production-expansion frozen roster differs: {model_kind}")
            latent[full_right] = np.asarray(v9_latent[v9_left], dtype=np.float32)

        frozen = cells["strict_reference"].to_numpy(dtype=bool) | cells["primary_query"].to_numpy(dtype=bool)
        before_frozen = np.asarray(latent[frozen]).copy()
        offsets = []
        cells["production_mapping_offset_source"] = "none"
        for dataset in secondary_datasets:
            result = _select_stage_blind_offset(latent, cells, dataset, minimum)
            if result is None:
                offsets.append({
                    "dataset": dataset, "applied": False,
                    "reason": "fewer_than_three_donors_in_a_required_fit_group",
                })
                continue
            apply_mask = result.pop("apply_mask")
            fit_mask = result.pop("fit_mask")
            before = np.asarray(latent[apply_mask], dtype=np.float64)
            latent[apply_mask] = (before + result["offset"]).astype(np.float32)
            after = np.asarray(latent[apply_mask], dtype=np.float64)
            centered_error = float(np.max(np.abs(
                (after - after.mean(axis=0)) - (before - before.mean(axis=0))
            ), initial=0.0))
            cells.loc[apply_mask, "production_mapping_offset_source"] = result["fit_scope"]
            result["offset"] = list(map(float, result["offset"]))
            result.update({
                "applied": True,
                "n_applied_cells": int(apply_mask.sum()),
                "n_fit_cells": int(fit_mask.sum()),
                "maximum_centered_coordinate_error": centered_error,
                "conditions_available_to_fit": False,
            })
            offsets.append(result)
        latent.flush()
        if not np.array_equal(np.asarray(latent[frozen]), before_frozen):
            raise ContractError(f"production expansion changed frozen V32 cells: {model_kind}")

        cells_path = model_output / "embedding_cells.tsv.gz"
        with DeterministicGzipTextWriter(cells_path) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(cells.columns)
            writer.writerows(cells.itertuples(index=False, name=None))
        embedding = {
            "schema_version": "masld-cl-embedding-v1",
            "config_sha256": config["_config_sha256"],
            "method": "replay_ewc_routed_stage_blind_production_expansion",
            "model_kind": model_kind,
            "n_cells": len(cells), "n_latent": 30,
            "n_analysis_ineligible": spec["projection_only_cells"],
            "latent_file": latent_path.name, "cells_file": cells_path.name,
            "latent_sha256": sha256_path(latent_path),
            "cells_sha256": sha256_path(cells_path),
        }
        write_json_exclusive(model_output / "embedding_manifest.json", embedding)
        manifest = {
            "schema_version": "masld-cl-production-expansion-model-v35",
            "config_sha256": config["_config_sha256"],
            "model_kind": model_kind,
            "lineage_label": spec["lineage_label"],
            "frozen_primary_reference_coordinates_bitwise_identical": True,
            "conditions_available_to_fit": False,
            "descriptive_only_cells_available_to_fit": False,
            "offsets": offsets,
            "embedding": embedding,
        }
        write_json_exclusive(model_output / "production_manifest.json", manifest)
        bundle_models[model_kind] = {
            "embedding_manifest": str((model_output / "embedding_manifest.json").resolve()),
            "embedding_manifest_sha256": sha256_path(model_output / "embedding_manifest.json"),
            "production_manifest": str((model_output / "production_manifest.json").resolve()),
            "production_manifest_sha256": sha256_path(model_output / "production_manifest.json"),
        }
    bundle = {
        "schema_version": "masld-cl-production-expansion-bundle-v35",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "secondary_datasets": secondary_datasets,
        "model_roster": list(policy["models"]),
        "models": bundle_models,
        "conditions_available_to_fit": False,
        "secondary_stress_status": "method_development_not_independent_confirmation",
    }
    write_json_exclusive(output_root / "bundle_manifest.json", bundle)
    return bundle
