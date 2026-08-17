"""Join locked raw-PCA geometry to replay-plus-EWC reference-only labels."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .embedding import load_embedding, matched_rows
from .firewall import validate_program_firewall
from .latent_knn_label_audit import _fit_predict_knn, load_latent_knn_label_audit_policy
from .training import _capped_indices


def build_routed_raw_geometry(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path = Path(policy_value).resolve()
    reference_root = Path(config["_config_path"]).resolve().parent / "reference"
    profiles = {
        reference_root / "replay_ewc_routed_raw_geometry_policy_v31.json": {
            "policy_schema": "masld-cl-replay-ewc-routed-raw-geometry-policy-v31",
            "method": "replay_ewc_routed_raw_geometry_adapter",
            "manifest_schema": "masld-cl-replay-ewc-routed-raw-geometry-v31",
            "geometry_name": "v14",
        },
        reference_root / "replay_ewc_routed_v9_geometry_policy_v32.json": {
            "policy_schema": "masld-cl-replay-ewc-routed-v9-geometry-policy-v32",
            "method": "replay_ewc_routed_v9_geometry_adapter",
            "manifest_schema": "masld-cl-replay-ewc-routed-v9-geometry-v32",
            "geometry_name": "v9",
        },
    }
    profile = profiles.get(policy_path)
    if profile is None:
        raise ContractError("routed geometry requires a source-controlled policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != profile["policy_schema"]
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("router", {}).get("query_labels_available_to_fit") is not False
        or policy.get("firewall", {}).get("no_new_hyperparameter_or_control_selection") is not True
    ):
        raise ContractError("routed-geometry policy identity differs")
    root = policy_path.parent
    sources = (
        ("geometry", "embedding", "embedding_sha256"),
        ("geometry", "selection_lock", "selection_lock_sha256"),
        ("geometry", "outcomes", "outcomes_sha256"),
        ("router", "embedding", "embedding_sha256"),
        ("router", "classifier_policy", "classifier_policy_sha256"),
        ("router", "v26_sensitivity", "v26_sensitivity_sha256"),
    )
    resolved = {}
    for section, key, hash_key in sources:
        path = (root / policy[section][key]).resolve()
        if sha256_path(path) != policy[section][hash_key]:
            raise ContractError(f"V31 policy source changed: {section}.{key}")
        resolved[(section, key)] = path
    if "seed_confirmation" in policy["geometry"]:
        path = (root / policy["geometry"]["seed_confirmation"]).resolve()
        if sha256_path(path) != policy["geometry"]["seed_confirmation_sha256"]:
            raise ContractError("V32 seed confirmation changed")
    with resolved[("geometry", "outcomes")].open() as handle:
        geometry_outcomes = json.load(handle)
    disease_pass = (
        geometry_outcomes.get("disease_preservation", {}).get("pass")
        if profile["geometry_name"] == "v14"
        else geometry_outcomes.get("disease_preservation_pass")
    )
    if disease_pass is not True:
        raise ContractError("routed source geometry did not pass disease preservation")
    with resolved[("geometry", "selection_lock")].open() as handle:
        geometry_selection = json.load(handle)
    if profile["geometry_name"] == "v9" and (
        float(geometry_selection.get("selected_global_reference_weight", -1))
        != policy["geometry"]["required_global_reference_weight"]
    ):
        raise ContractError("V32 source geometry is not the locked V9 weight")

    geometry_info, geometry_latent, cells = load_embedding(
        resolved[("geometry", "embedding")]
    )
    _, router_latent, router_cells = load_embedding(resolved[("router", "embedding")])
    router_left, geometry_right = matched_rows(router_cells, cells)
    order = np.argsort(geometry_right)
    router_left, geometry_right = router_left[order], geometry_right[order]
    if not np.array_equal(geometry_right, np.arange(len(cells))):
        raise ContractError("V31 geometry and router rosters differ")
    ordered_router_latent = np.asarray(router_latent)[router_left]
    ordered_router_cells = router_cells.iloc[router_left].reset_index(drop=True)
    if not np.array_equal(
        ordered_router_cells["cell_id"].astype(str).to_numpy(),
        cells["cell_id"].astype(str).to_numpy(),
    ):
        raise ContractError("V31 router cell order differs")
    _, classifier_policy = load_latent_knn_label_audit_policy(
        config, resolved[("router", "classifier_policy")]
    )
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    query = cells["analysis_eligible"].to_numpy(dtype=bool) & ~reference
    reference_pool = _capped_indices(
        type("ADataView", (), {"obs": cells, "n_obs": len(cells)})(),
        np.flatnonzero(reference), "all_lineage", config,
        classifier_policy["classifier"]["reference_cap_seed"],
    )
    predictions, routing_identity = _fit_predict_knn(
        ordered_router_latent[reference_pool],
        cells.iloc[reference_pool]["audit_cell_type"].astype(str).to_numpy(),
        ordered_router_latent[query], classifier_policy["classifier"],
    )
    output_cells = cells.copy()
    output_cells["routing_label"] = output_cells["predicted_cell_type"].astype(str)
    output_cells.loc[query, "routing_label"] = predictions

    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    source_latent_path = resolved[("geometry", "embedding")].parent / geometry_info["latent_file"]
    latent_path = output / geometry_info["latent_file"]
    latent_path.hardlink_to(source_latent_path)
    if sha256_path(latent_path) != geometry_info["latent_sha256"]:
        raise ContractError("V31 geometry hard link differs")
    cells_path = output / geometry_info["cells_file"]
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_cells.columns)
        writer.writerows(output_cells.itertuples(index=False, name=None))
    embedding = dict(geometry_info)
    embedding.update({
        "method": profile["method"],
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
    })
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    manifest = {
        "schema_version": profile["manifest_schema"],
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "geometry_bitwise_identical_to_source": True,
        "geometry_source": profile["geometry_name"],
        "query_labels_available_to_fit": False,
        "routing_identity": routing_identity,
        "embedding": embedding,
    }
    write_json_exclusive(output / "adapter_manifest.json", manifest)
    return manifest
