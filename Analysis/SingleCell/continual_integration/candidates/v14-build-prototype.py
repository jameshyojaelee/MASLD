#!/usr/bin/env python3
"""Build the prospectively locked V14 all-lineage control-only prototype grid."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from masld_cl.config import load_config, write_json_exclusive
from masld_cl.contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from masld_cl.control_adapter import _donor_balanced_centroid
from masld_cl.embedding import load_embedding, matched_rows
from masld_cl.orthogonal_bridge import (
    _balanced_fit_positions,
    apply_scaled_orthogonal_bridge,
    fit_scaled_orthogonal_bridge,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--raw-pca-embedding", required=True)
    parser.add_argument("--parent-embedding", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    policy_path = Path(args.policy).resolve()
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version")
        != "masld-cl-geometry-preserving-continual-adapter-policy-v14"
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("firewall", {}).get(
            "no_query_label_disease_or_held_study_metric_may_be_read_before_the_v14_control_lock"
        )
        is not True
    ):
        raise ContractError("V14 policy identity or firewall differs")

    root = Path(config["_config_path"]).resolve().parent
    parent_path = Path(args.parent_embedding).resolve()
    expected_parent = root / policy["v13_parent"]["embedding_relative_path"].removeprefix(
        "candidates/"
    )
    expected_parent = root / "candidates" / Path(
        policy["v13_parent"]["embedding_relative_path"]
    ).relative_to("candidates")
    if (
        parent_path != expected_parent
        or sha256_path(parent_path) != policy["v13_parent"]["embedding_sha256"]
    ):
        raise ContractError("V14 parent differs from the prospective policy")
    parent = load_embedding(parent_path)
    parent_run_path = parent_path.parent / "update_manifest.json"
    with parent_run_path.open() as handle:
        parent_run = json.load(handle)
    if (
        parent_run.get("method") != "continual_distillation"
        or parent_run.get("ewc_lambda") != 100.0
        or parent_run.get("replay_fraction") != 0.2
        or parent_run.get("distillation_weight") != 100.0
        or parent_run.get("embedding") != parent[0]
    ):
        raise ContractError("V14 parent training contract differs")

    raw_path = Path(args.raw_pca_embedding).resolve()
    reference_path = Path(args.reference_embedding).resolve()
    raw = load_embedding(raw_path)
    reference = load_embedding(reference_path)
    if raw[0].get("method") != "uncorrected_scalesc_pca_30d":
        raise ContractError("V14 raw-PCA source differs")
    raw_left, parent_right = matched_rows(raw[2], parent[2])
    order = np.argsort(parent_right)
    raw_left, parent_right = raw_left[order], parent_right[order]
    if not np.array_equal(parent_right, np.arange(len(parent[2]))):
        raise ContractError("raw PCA and V13 parent cell rosters differ")
    raw_values = np.asarray(raw[1])[raw_left]
    cells = parent[2].copy()

    reference_positions = np.flatnonzero(cells["strict_reference"].to_numpy(dtype=bool))
    ref_left, parent_ref_right = matched_rows(
        reference[2], cells.iloc[reference_positions].reset_index(drop=True)
    )
    parent_reference_positions = reference_positions[parent_ref_right]
    fit_positions = _balanced_fit_positions(
        cells, parent_reference_positions, "all_lineage", config
    )
    ref_by_cell = {
        cell: index for index, cell in enumerate(reference[2]["cell_id"].astype(str))
    }
    fit_reference_positions = np.asarray(
        [ref_by_cell[cell] for cell in cells.iloc[fit_positions]["cell_id"].astype(str)],
        dtype=np.int64,
    )
    source_mean, target_mean, rotation, scale = fit_scaled_orthogonal_bridge(
        raw_values[fit_positions], np.asarray(reference[1])[fit_reference_positions]
    )
    bridged = apply_scaled_orthogonal_bridge(
        raw_values, source_mean, target_mean, rotation, scale
    )
    bridged[parent_reference_positions] = np.asarray(reference[1])[ref_left]

    is_reference = cells["strict_reference"].to_numpy(dtype=bool)
    is_primary = cells["primary_query"].to_numpy(dtype=bool)
    is_control = is_primary & cells["query_control"].to_numpy(dtype=bool) & ~is_reference
    preparations = cells["preparation"].astype(str).to_numpy()
    datasets = cells["dataset"].astype(str).to_numpy()
    global_centroid, n_global = _donor_balanced_centroid(
        bridged, cells, is_reference
    )
    compatible = is_reference & (
        preparations == policy["control_calibration"]["compatible_reference_preparation"]
    )
    compatible_centroid, n_compatible = _donor_balanced_centroid(
        bridged, cells, compatible
    )
    if n_compatible < 3:
        raise ContractError("V14 compatible reference group is underpowered")

    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=False)
    outputs = []
    for global_weight in map(
        float, policy["control_calibration"]["global_reference_weight_grid"]
    ):
        target = (
            (1.0 - global_weight) * compatible_centroid
            + global_weight * global_centroid
        )
        for offset_weight in map(
            float, policy["control_calibration"]["control_offset_weight_grid"]
        ):
            adapted = np.asarray(bridged).copy()
            offsets = []
            for dataset in policy["control_calibration"]["primary_query_datasets"]:
                control_mask = is_control & (datasets == dataset)
                control_centroid, n_control = _donor_balanced_centroid(
                    bridged, cells, control_mask
                )
                if n_control < int(
                    policy["control_calibration"]["minimum_control_donors_per_dataset"]
                ):
                    raise ContractError(f"underpowered V14 control group: {dataset}")
                offset = offset_weight * (target - control_centroid)
                apply_mask = is_primary & (datasets == dataset)
                before = np.asarray(bridged[apply_mask], dtype=np.float64)
                adapted[apply_mask] = (before + offset).astype(adapted.dtype, copy=False)
                after = np.asarray(adapted[apply_mask], dtype=np.float64)
                centered_error = float(
                    np.max(
                        np.abs(
                            (after - after.mean(axis=0))
                            - (before - before.mean(axis=0))
                        ),
                        initial=0.0,
                    )
                )
                offsets.append(
                    {
                        "dataset": dataset,
                        "n_control_donors": n_control,
                        "offset": list(map(float, offset)),
                        "maximum_centered_coordinate_error": centered_error,
                    }
                )
            if not np.array_equal(adapted[is_reference], bridged[is_reference]):
                raise ContractError("V14 calibration changed reference coordinates")

            setting = f"g{int(round(global_weight * 100)):03d}_w{int(round(offset_weight * 100)):03d}"
            output = output_root / setting
            output.mkdir()
            latent_path = output / "embedding_latent.npy"
            cells_path = output / parent[0]["cells_file"]
            np.save(latent_path, adapted)
            with DeterministicGzipTextWriter(cells_path) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(cells.columns)
                writer.writerows(cells.itertuples(index=False, name=None))
            embedding = dict(parent[0])
            embedding.update(
                {
                    "method": "geometry_preserving_continual_adapter",
                    "latent_file": latent_path.name,
                    "cells_file": cells_path.name,
                    "latent_sha256": sha256_path(latent_path),
                    "cells_sha256": sha256_path(cells_path),
                }
            )
            write_json_exclusive(output / "embedding_manifest.json", embedding)
            manifest = {
                "schema_version": "masld-cl-geometry-preserving-continual-adapter-v14",
                "config_sha256": config["_config_sha256"],
                "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
                "control_only": True,
                "outcomes_unlocked": False,
                "model_kind": "all_lineage",
                "global_reference_weight": global_weight,
                "control_offset_weight": offset_weight,
                "n_global_reference_donors": n_global,
                "n_compatible_reference_donors": n_compatible,
                "reference_coordinates_bitwise_frozen": True,
                "query_predictions_source": str(parent_path),
                "query_predictions_source_sha256": sha256_path(parent_path),
                "raw_pca_source": str(raw_path),
                "raw_pca_source_sha256": sha256_path(raw_path),
                "reference_source": str(reference_path),
                "reference_source_sha256": sha256_path(reference_path),
                "bridge": {
                    "fit_cells": int(len(fit_positions)),
                    "scale": scale,
                    "source_mean": source_mean.tolist(),
                    "target_mean": target_mean.tolist(),
                    "rotation": rotation.tolist(),
                },
                "offsets": offsets,
                "embedding": embedding,
            }
            write_json_exclusive(output / "adapter_manifest.json", manifest)
            outputs.append(
                {
                    "setting": setting,
                    "global_reference_weight": global_weight,
                    "control_offset_weight": offset_weight,
                    "embedding": str((output / "embedding_manifest.json").resolve()),
                    "manifest": str((output / "adapter_manifest.json").resolve()),
                }
            )
    grid = {
        "schema_version": "masld-cl-geometry-preserving-continual-adapter-grid-v14",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "control_only": True,
        "outcomes_unlocked": False,
        "outputs": outputs,
    }
    write_json_exclusive(output_root / "grid_manifest.json", grid)
    print(json.dumps(grid, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
