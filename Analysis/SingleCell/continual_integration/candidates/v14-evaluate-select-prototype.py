#!/usr/bin/env python3
"""Evaluate controls and freeze the prospective V14 all-lineage pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from masld_cl.backbone_diagnostic import _bundle_centroids, _paired_improvement
from masld_cl.config import canonical_json_bytes, load_config, write_json_exclusive
from masld_cl.contracts import ContractError, sha256_path
from masld_cl.control_evaluation import (
    _neighbor_jaccard_loss,
    _reference_f1_change_ci,
    _shift_and_standard_error,
)
from masld_cl.embedding import load_embedding, matched_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--grid", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--architecture-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    policy_path = Path(args.policy).resolve()
    grid_path = Path(args.grid).resolve()
    with policy_path.open() as handle:
        policy = json.load(handle)
    with grid_path.open() as handle:
        grid = json.load(handle)
    if (
        policy.get("schema_version")
        != "masld-cl-geometry-preserving-continual-adapter-policy-v14"
        or policy.get("config_sha256") != config["_config_sha256"]
        or grid.get("schema_version")
        != "masld-cl-geometry-preserving-continual-adapter-grid-v14"
        or grid.get("policy", {}).get("sha256") != sha256_path(policy_path)
        or grid.get("outcomes_unlocked") is not False
        or len(grid.get("outputs", [])) != 15
    ):
        raise ContractError("V14 policy or grid identity differs")

    reference_path = Path(args.reference_embedding).resolve()
    harmony_path = Path(args.harmony_embedding).resolve()
    architecture_path = Path(args.architecture_embedding).resolve()
    reference = load_embedding(reference_path)
    harmony = load_embedding(harmony_path)
    architecture = load_embedding(architecture_path)
    if harmony[0].get("method") != "incumbent_scalesc_pca_harmony":
        raise ContractError("V14 Harmony comparator differs")

    results = []
    for index, source in enumerate(grid["outputs"]):
        embedding_path = Path(source["embedding"]).resolve()
        manifest_path = Path(source["manifest"]).resolve()
        candidate = load_embedding(embedding_path)
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        if (
            manifest.get("schema_version")
            != "masld-cl-geometry-preserving-continual-adapter-v14"
            or manifest.get("policy", {}).get("sha256") != sha256_path(policy_path)
            or manifest.get("outcomes_unlocked") is not False
            or manifest.get("embedding") != candidate[0]
        ):
            raise ContractError(f"V14 candidate manifest differs: {source['setting']}")

        positions = np.flatnonzero(candidate[2]["strict_reference"].to_numpy(dtype=bool))
        candidate_reference = candidate[2].iloc[positions].reset_index(drop=True)
        left, right = matched_rows(reference[2], candidate_reference)
        reference_cells = reference[2].iloc[left].reset_index(drop=True)
        candidate_cells = candidate_reference.iloc[right].reset_index(drop=True)
        reference_latent = np.asarray(reference[1])[left]
        candidate_latent = np.asarray(candidate[1])[positions[right]]
        exact_coordinates = np.array_equal(reference_latent, candidate_latent)
        exact_predictions = np.array_equal(
            reference_cells["predicted_cell_type"].astype(str).to_numpy(),
            candidate_cells["predicted_cell_type"].astype(str).to_numpy(),
        )
        if not exact_coordinates or not exact_predictions:
            raise ContractError(f"V14 reference export is not exact: {source['setting']}")
        f1 = _reference_f1_change_ci(
            reference_cells,
            candidate_cells,
            config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 1901 + index,
        )
        jaccard = _neighbor_jaccard_loss(
            reference_latent,
            candidate_latent,
            reference_cells,
            k=config["evaluation"]["reference_neighborhood_k"],
            maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
            seed=config["screen"]["seed"] + 1901 + index,
        )
        harmony_result = _paired_improvement(
            _bundle_centroids(candidate),
            _bundle_centroids(harmony),
            config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 2001 + index,
        )
        architecture_result = _paired_improvement(
            _bundle_centroids(candidate),
            _bundle_centroids(architecture),
            config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 2101 + index,
        )
        shift, shift_se = _shift_and_standard_error(
            candidate[2],
            candidate[1],
            config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 2201 + index,
        )
        gates = {
            "reference_macro_f1": f1["ci_low"]
            > policy["selection"]["required_reference_macro_f1_change_ci_low"],
            "reference_neighborhood": jaccard
            <= policy["selection"]["maximum_reference_neighborhood_jaccard_loss"],
            "alignment_vs_harmony": harmony_result["improvement"]
            >= policy["selection"]["minimum_alignment_improvement_vs_harmony"]
            and harmony_result["ci_low"]
            > policy["selection"]["alignment_ci_low_must_exceed"],
            "alignment_vs_architecture_surgery": architecture_result["improvement"]
            >= config["gates"]["minimum_control_alignment_improvement"]
            and architecture_result["ci_low"] > 0,
        }
        results.append(
            {
                "setting": source["setting"],
                "global_reference_weight": source["global_reference_weight"],
                "control_offset_weight": source["control_offset_weight"],
                "reference_macro_f1_change": f1,
                "reference_neighborhood_jaccard_loss": jaccard,
                "reference_coordinates_bitwise_exact": exact_coordinates,
                "reference_predictions_bitwise_exact": exact_predictions,
                "shift_control": shift,
                "shift_control_standard_error": shift_se,
                "alignment_vs_harmony": harmony_result,
                "alignment_vs_architecture_surgery": architecture_result,
                "gates": gates,
                "eligible": all(gates.values()),
                "embedding": str(embedding_path),
                "embedding_sha256": sha256_path(embedding_path),
                "manifest": str(manifest_path),
                "manifest_sha256": sha256_path(manifest_path),
            }
        )

    survivors = [row for row in results if row["eligible"]]
    if not survivors:
        selected = None
        threshold = None
    else:
        best = min(survivors, key=lambda row: row["shift_control"])
        threshold = best["shift_control"] + best["shift_control_standard_error"]
        within_one_se = [
            row for row in survivors if row["shift_control"] <= threshold
        ]
        within_one_se.sort(
            key=lambda row: (
                row["control_offset_weight"],
                abs(row["global_reference_weight"] - 0.75),
                row["global_reference_weight"],
                row["shift_control"],
            )
        )
        selected = within_one_se[0]

    decision = {
        "schema_version": "masld-cl-geometry-preserving-continual-adapter-selection-v14",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "grid": {"path": str(grid_path), "sha256": sha256_path(grid_path)},
        "control_only": True,
        "outcome_variables_read": [],
        "selection_frozen": selected is not None,
        "outcomes_unlocked": selected is not None,
        "selected": selected,
        "one_standard_error_threshold": threshold,
        "eligible_settings": [row["setting"] for row in survivors],
        "results": results,
        "comparators": {
            "reference": {"path": str(reference_path), "sha256": sha256_path(reference_path)},
            "harmony": {"path": str(harmony_path), "sha256": sha256_path(harmony_path)},
            "architecture_surgery": {
                "path": str(architecture_path),
                "sha256": sha256_path(architecture_path),
            },
        },
    }
    decision["lock_sha256"] = hashlib.sha256(canonical_json_bytes(decision)).hexdigest()
    write_json_exclusive(args.output, decision)
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
