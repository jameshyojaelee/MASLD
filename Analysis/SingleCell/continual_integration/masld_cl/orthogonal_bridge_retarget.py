"""Control-only reference-target retuning of the geometry-preserving bridge."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .control_adapter import _donor_balanced_centroid
from .embedding import load_embedding
from .execution import source_identity, verify_source_identity_payload
from .firewall import validate_program_firewall


def build_retargeted_bridge_grid(
    config: dict[str, Any], policy_value: str | Path,
    parent_embedding_value: str | Path, output_root_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    policy_path = Path(policy_value).resolve()
    if policy_path != pipeline_root / "reference" / "orthogonal_bridge_retarget_policy_v9.json":
        raise ContractError("bridge retargeting requires the source-controlled V9 policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != "masld-cl-orthogonal-bridge-retarget-policy-v9"
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("case_stage_program_hero_gene_and_cas13_outcomes_locked") is not True
    ):
        raise ContractError("bridge retarget policy is invalid")
    parent_results = []
    for source in policy["parent_control_results"]:
        path = pipeline_root / source["relative_path"]
        if sha256_path(path) != source["sha256"]:
            raise ContractError("bridge retarget parent control result changed")
        with path.open() as handle:
            result = json.load(handle)
        if result.get("outcomes_unlocked") is not False:
            raise ContractError("bridge retarget parent exposed outcomes")
        parent_results.append(result)
    if {row["model_kind"] for row in parent_results} != {"all_lineage", *config["lineages"]}:
        raise ContractError("bridge retarget parent model roster is incomplete")

    parent_path = Path(parent_embedding_value).resolve()
    info, latent, cells = load_embedding(parent_path)
    parent_manifest_path = parent_path.parent / "orthogonal_bridge_manifest.json"
    with parent_manifest_path.open() as handle:
        parent = json.load(handle)
    verify_source_identity_payload(parent.get("source_identity"))
    if (
        parent.get("schema_version") != "masld-cl-orthogonal-bridge-v8"
        or parent.get("embedding") != info
        or float(parent.get("control_offset_weight", -1)) != float(policy["parent_control_offset_weight"])
        or parent.get("model_kind") not in {"all_lineage", *config["lineages"]}
    ):
        raise ContractError("bridge retarget parent candidate is invalid")
    matching = [row for row in parent_results if row["model_kind"] == parent["model_kind"]]
    if len(matching) != 1 or Path(matching[0]["sources"]["candidate_embedding"]).resolve() != parent_path:
        raise ContractError("bridge retarget parent was not in the locked control matrix")

    reference = cells["strict_reference"].to_numpy(dtype=bool)
    primary = cells["primary_query"].to_numpy(dtype=bool)
    controls = primary & cells["query_control"].to_numpy(dtype=bool) & ~reference
    preparations = cells["preparation"].astype(str).to_numpy()
    datasets = cells["dataset"].astype(str).to_numpy()
    global_centroid, n_global = _donor_balanced_centroid(latent, cells, reference)
    compatible = reference & (preparations == policy["compatible_reference_preparation"])
    compatible_centroid, n_compatible = _donor_balanced_centroid(latent, cells, compatible)
    if n_compatible < 3:
        raise ContractError("bridge retarget compatible reference is underpowered")

    output_root = Path(output_root_value)
    output_root.mkdir(parents=True, exist_ok=False)
    frozen_source = source_identity(pipeline_root)
    outputs = []
    for weight in map(float, policy["global_reference_weight_grid"]):
        target = (1.0 - weight) * compatible_centroid + weight * global_centroid
        adapted = np.asarray(latent).copy()
        offsets = []
        for dataset in policy["primary_query_datasets"]:
            control_mask = controls & (datasets == dataset)
            control_centroid, n_control = _donor_balanced_centroid(latent, cells, control_mask)
            offset = target - control_centroid
            apply_mask = primary & (datasets == dataset)
            before = np.asarray(latent[apply_mask], dtype=np.float64)
            adapted[apply_mask] = (before + offset).astype(adapted.dtype, copy=False)
            after = np.asarray(adapted[apply_mask], dtype=np.float64)
            centered_error = float(np.max(np.abs(
                (after - after.mean(axis=0)) - (before - before.mean(axis=0))
            ), initial=0.0))
            offsets.append({
                "dataset": dataset, "n_control_donors": n_control,
                "n_global_reference_donors": n_global,
                "n_compatible_reference_donors": n_compatible,
                "global_reference_weight": weight, "offset": list(map(float, offset)),
                "maximum_centered_coordinate_error": centered_error,
            })
        if not np.array_equal(adapted[reference], np.asarray(latent)[reference]):
            raise ContractError("bridge retargeting changed frozen reference coordinates")
        output = output_root / f"g{int(round(weight * 100)):03d}"
        output.mkdir()
        latent_path = output / "embedding_latent.npy"
        cells_path = output / info["cells_file"]
        np.save(latent_path, adapted)
        with DeterministicGzipTextWriter(cells_path) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(cells.columns)
            writer.writerows(cells.itertuples(index=False, name=None))
        embedding = dict(info)
        embedding.update({
            "method": "retargeted_orthogonal_bridge_continual_adapter",
            "latent_file": latent_path.name, "cells_file": cells_path.name,
            "latent_sha256": sha256_path(latent_path), "cells_sha256": sha256_path(cells_path),
        })
        write_json_exclusive(output / "embedding_manifest.json", embedding)
        manifest = {
            "schema_version": "masld-cl-orthogonal-bridge-retarget-v9",
            "config_sha256": config["_config_sha256"], "policy_realpath": str(policy_path),
            "policy_sha256": sha256_path(policy_path), "source_identity": frozen_source,
            "model_kind": parent["model_kind"], "control_offset_weight": 1.0,
            "global_reference_weight": weight, "control_only": True,
            "outcomes_unlocked": False,
            "reference_coordinates_and_predictions_frozen_to_preupdate": True,
            "query_predictions_frozen_to_architecture_surgery": True,
            "parent_embedding_realpath": str(parent_path),
            "parent_embedding_sha256": sha256_path(parent_path),
            "parent_manifest_realpath": str(parent_manifest_path),
            "parent_manifest_sha256": sha256_path(parent_manifest_path),
            "reference_embedding_realpath": parent["reference_embedding_realpath"],
            "architecture_embedding_realpath": parent["architecture_embedding_realpath"],
            "offsets": offsets, "embedding": embedding,
        }
        write_json_exclusive(output / "orthogonal_bridge_manifest.json", manifest)
        outputs.append({
            "global_reference_weight": weight,
            "embedding_manifest": str((output / "embedding_manifest.json").resolve()),
            "bridge_manifest": str((output / "orthogonal_bridge_manifest.json").resolve()),
            "bridge_manifest_sha256": sha256_path(output / "orthogonal_bridge_manifest.json"),
        })
    grid = {
        "schema_version": "masld-cl-orthogonal-bridge-retarget-grid-v9",
        "config_sha256": config["_config_sha256"], "model_kind": parent["model_kind"],
        "policy_sha256": sha256_path(policy_path), "source_identity": frozen_source,
        "outputs": outputs,
    }
    write_json_exclusive(output_root / "grid_manifest.json", grid)
    return grid
