"""Project the confirmed V32 candidate over the immutable full atlas roster."""

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


POLICY_SCHEMA = "masld-cl-full-atlas-projection-policy-v33"
MANIFEST_SCHEMA = "masld-cl-full-atlas-projection-v33"


def load_full_atlas_projection_policy(
    config: dict[str, Any], value: str | Path,
) -> tuple[Path, dict[str, Any], dict[str, Path]]:
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "full_atlas_projection_policy_v33.json"
    )
    if path != expected:
        raise ContractError("full-atlas projection requires the source-controlled V33 policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("projection", {}).get("secondary_query_labels_available_to_fit") is not False
        or policy.get("firewall", {}).get("no_refitting_or_hyperparameter_selection") is not True
        or policy.get("firewall", {}).get("secondary_and_descriptive_cells_do_not_steer_parameters") is not True
    ):
        raise ContractError("full-atlas projection policy identity differs")
    resolved: dict[str, Path] = {}
    for name, source in policy.get("sources", {}).items():
        source_path = (path.parent / source["path"]).resolve()
        if sha256_path(source_path) != source["sha256"]:
            raise ContractError(f"full-atlas projection source changed: {name}")
        resolved[name] = source_path
    required = {
        "full_raw_pca", "v8_bridge", "v9_retarget", "v32_embedding",
        "v32_confirmation", "classifier_policy",
    }
    if set(resolved) != required:
        raise ContractError("full-atlas projection source roster differs")
    with resolved["v32_confirmation"].open() as handle:
        confirmation = json.load(handle)
    if confirmation.get("confirmation_pass") is not True:
        raise ContractError("V32 confirmation is not passing")
    return path, policy, resolved


def _project_frozen_bridge(
    raw: np.ndarray, bridge: dict[str, Any], output: np.ndarray, chunk_size: int = 100_000,
) -> None:
    parameters = bridge.get("bridge", {})
    source_mean = np.asarray(parameters.get("source_mean"), dtype=np.float64)
    target_mean = np.asarray(parameters.get("target_mean"), dtype=np.float64)
    rotation = np.asarray(parameters.get("rotation"), dtype=np.float64)
    scale = float(parameters.get("scale", np.nan))
    dimensions = raw.shape[1]
    if (
        source_mean.shape != (dimensions,)
        or target_mean.shape != (dimensions,)
        or rotation.shape != (dimensions, dimensions)
        or not np.isfinite(scale)
        or scale <= 0
        or np.max(np.abs(rotation.T @ rotation - np.eye(dimensions))) > 1e-10
    ):
        raise ContractError("frozen bridge parameters are invalid")
    for start in range(0, len(raw), chunk_size):
        stop = min(start + chunk_size, len(raw))
        values = np.asarray(raw[start:stop], dtype=np.float64)
        output[start:stop] = (
            (values - source_mean) @ rotation * scale + target_mean
        ).astype(np.float32)


def build_full_atlas_projection(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy, sources = load_full_atlas_projection_policy(config, policy_value)
    expected = policy["expected_contract"]

    raw_info, raw_latent, cells = load_embedding(sources["full_raw_pca"])
    v32_info, v32_latent, v32_cells = load_embedding(sources["v32_embedding"])
    counts = {
        "descriptive_cells": len(cells),
        "analyzed_cells": int(cells["analysis_eligible"].sum()),
        "strict_reference_cells": int(cells["strict_reference"].sum()),
        "primary_query_cells": int(cells["primary_query"].sum()),
        "secondary_mapping_cells": int((
            cells["analysis_eligible"] & ~cells["strict_reference"] & ~cells["primary_query"]
        ).sum()),
        "descriptive_only_cells": int((~cells["analysis_eligible"]).sum()),
        "v32_cells": len(v32_cells),
        "latent_dimensions": int(raw_latent.shape[1]),
    }
    if counts != expected:
        raise ContractError(f"full-atlas projection contract differs: {counts}")
    if raw_info.get("method") != "uncorrected_scalesc_pca_30d":
        raise ContractError("V33 raw-PCA source method differs")

    with sources["v8_bridge"].open() as handle:
        bridge = json.load(handle)
    with sources["v9_retarget"].open() as handle:
        retarget = json.load(handle)
    if (
        bridge.get("schema_version") != "masld-cl-orthogonal-bridge-v8"
        or bridge.get("raw_pca_embedding_sha256") != sha256_path(sources["full_raw_pca"])
        or retarget.get("schema_version") != "masld-cl-orthogonal-bridge-retarget-v9"
        or float(retarget.get("global_reference_weight", -1)) != 0.75
        or retarget.get("parent_manifest_sha256") != sha256_path(sources["v8_bridge"])
    ):
        raise ContractError("V33 frozen geometry chain differs")

    v32_left, full_right = matched_rows(v32_cells, cells)
    order = np.argsort(full_right)
    v32_left, full_right = v32_left[order], full_right[order]
    if len(full_right) != expected["v32_cells"]:
        raise ContractError("V33 does not recover the complete V32 roster")
    if not np.array_equal(
        v32_cells.iloc[v32_left]["cell_id"].astype(str).to_numpy(),
        cells.iloc[full_right]["cell_id"].astype(str).to_numpy(),
    ):
        raise ContractError("V33 V32 cell matching differs")

    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    latent_path = output / "embedding_latent.npy"
    projected = np.lib.format.open_memmap(
        latent_path, mode="w+", dtype=np.float32,
        shape=(expected["descriptive_cells"], expected["latent_dimensions"]),
    )
    _project_frozen_bridge(raw_latent, bridge, projected)
    # The selected primary/reference candidate is the frozen authority. Copying
    # these coordinates also avoids platform-level BLAS drift in the projection.
    projected[full_right] = np.asarray(v32_latent[v32_left], dtype=np.float32)
    projected.flush()
    if not np.array_equal(projected[full_right], np.asarray(v32_latent[v32_left])):
        raise ContractError("V33 failed to reproduce V32 coordinates bitwise")

    output_cells = cells.copy()
    output_cells["routing_label"] = ""
    output_cells["routing_source"] = ""
    v32_labels = v32_cells.iloc[v32_left]["routing_label"].astype(str).to_numpy()
    output_cells.loc[full_right, "routing_label"] = v32_labels
    output_cells.loc[full_right, "routing_source"] = "v32_replay_ewc_reference_only"
    remaining = np.ones(len(cells), dtype=bool)
    remaining[full_right] = False
    remaining_positions = np.flatnonzero(remaining)
    if (
        int((~cells.iloc[remaining_positions]["analysis_eligible"]).sum())
        != expected["descriptive_only_cells"]
    ):
        raise ContractError("V33 descriptive roster is not confined to frozen projection cells")
    _, classifier_policy = load_latent_knn_label_audit_policy(
        config, sources["classifier_policy"]
    )
    reference_positions = np.flatnonzero(cells["strict_reference"].to_numpy(dtype=bool))
    reference_pool = _capped_indices(
        type("ADataView", (), {"obs": cells, "n_obs": len(cells)})(),
        reference_positions, "all_lineage", config,
        classifier_policy["classifier"]["reference_cap_seed"],
    )
    predictions, routing_identity = _fit_predict_knn(
        projected[reference_pool],
        cells.iloc[reference_pool]["audit_cell_type"].astype(str).to_numpy(),
        projected[remaining_positions], classifier_policy["classifier"],
    )
    output_cells.loc[remaining_positions, "routing_label"] = predictions
    output_cells.loc[remaining_positions, "routing_source"] = "v33_frozen_geometry_reference_only"
    if (output_cells["routing_label"].astype(str).str.len() == 0).any():
        raise ContractError("V33 produced an empty routing label")

    cells_path = output / "embedding_cells.tsv.gz"
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_cells.columns)
        writer.writerows(output_cells.itertuples(index=False, name=None))
    embedding = dict(raw_info)
    embedding.update({
        "method": "replay_ewc_routed_v9_full_atlas_projection",
        "latent_file": latent_path.name,
        "cells_file": cells_path.name,
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
        "n_analysis_ineligible": expected["descriptive_only_cells"],
    })
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    manifest = {
        "schema_version": MANIFEST_SCHEMA,
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "contract": counts,
        "v32_coordinates_bitwise_identical": True,
        "secondary_and_descriptive_cells_used_for_fit": False,
        "secondary_control_offsets_applied": False,
        "routing_identity": routing_identity,
        "routing_counts": {
            str(key): int(value) for key, value in
            output_cells["routing_source"].value_counts().sort_index().items()
        },
        "embedding": embedding,
    }
    write_json_exclusive(output / "projection_manifest.json", manifest)
    return manifest
