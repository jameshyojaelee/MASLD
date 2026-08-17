#!/usr/bin/env python3
"""Control-only evaluation for the locked V13 all-lineage pilot grid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from masld_cl.backbone_diagnostic import _bundle_centroids, _paired_improvement
from masld_cl.config import load_config, write_json_exclusive
from masld_cl.contracts import ContractError, sha256_path
from masld_cl.control_evaluation import (
    _neighbor_jaccard_loss,
    _reference_f1_change_ci,
)
from masld_cl.embedding import load_embedding, matched_rows


def _candidate_mapping(values: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        setting, separator, path = value.partition("=")
        if separator != "=" or not setting or setting in result:
            raise ContractError(f"invalid or duplicate candidate mapping: {value}")
        result[setting] = Path(path).resolve()
    if not result:
        raise ContractError("at least one candidate embedding is required")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--candidate", action="append", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    policy_path = Path(args.policy).resolve()
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != "masld-cl-continual-distillation-policy-v13"
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("selection_metrics") != [
            "reference_macro_f1_change_ci_low",
            "reference_neighborhood_jaccard_loss",
            "donor_balanced_query_control_shift",
            "alignment_improvement_vs_harmony",
        ]
    ):
        raise ContractError("V13 policy identity or selection metrics differ")

    reference_path = Path(args.reference_embedding).resolve()
    harmony_path = Path(args.harmony_embedding).resolve()
    reference = load_embedding(reference_path)
    harmony = load_embedding(harmony_path)
    if harmony[0].get("method") != "incumbent_scalesc_pca_harmony":
        raise ContractError("Harmony comparator role differs")

    results = {}
    for index, (setting, candidate_path) in enumerate(
        sorted(_candidate_mapping(args.candidate).items())
    ):
        candidate = load_embedding(candidate_path)
        run_path = candidate_path.parent / "update_manifest.json"
        with run_path.open() as handle:
            run = json.load(handle)
        if (
            run.get("method") != "continual_distillation"
            or run.get("seed") != 17
            or run.get("replay_fraction") != 0.2
            or run.get("distillation_weight") != 100.0
            or run.get("embedding", {}).get(
                "reference_coordinates_and_predictions_frozen"
            )
            is not True
        ):
            raise ContractError(f"candidate differs from the locked V13 method: {setting}")

        positions = np.flatnonzero(candidate[2]["strict_reference"].to_numpy(dtype=bool))
        candidate_reference = candidate[2].iloc[positions].reset_index(drop=True)
        left, right = matched_rows(reference[2], candidate_reference)
        if len(left) != len(reference[2]):
            raise ContractError(f"candidate lacks strict-reference cells: {setting}")
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
            raise ContractError(f"frozen reference export is not bitwise exact: {setting}")
        f1 = _reference_f1_change_ci(
            reference_cells,
            candidate_cells,
            config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 1701 + index,
        )
        jaccard = _neighbor_jaccard_loss(
            reference_latent,
            candidate_latent,
            reference_cells,
            k=config["evaluation"]["reference_neighborhood_k"],
            maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
            seed=config["screen"]["seed"] + 1701 + index,
        )
        alignment = _paired_improvement(
            _bundle_centroids(candidate),
            _bundle_centroids(harmony),
            config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 1801 + index,
        )
        gates = {
            "reference_macro_f1": f1["ci_low"]
            > -config["gates"]["reference_macro_f1_margin"],
            "reference_neighborhood": jaccard
            <= config["gates"]["reference_neighborhood_jaccard_loss"],
            "alignment_vs_harmony": alignment["improvement"]
            >= config["gates"]["minimum_control_alignment_improvement"]
            and alignment["ci_low"] > 0,
        }
        results[setting] = {
            "ewc_lambda": float(run["ewc_lambda"]),
            "stopping_epoch": int(run["stopping_epoch"]),
            "reference_macro_f1_change": f1,
            "reference_neighborhood_jaccard_loss": jaccard,
            "reference_coordinates_bitwise_exact": exact_coordinates,
            "reference_predictions_bitwise_exact": exact_predictions,
            "alignment_vs_harmony": alignment,
            "gates": gates,
            "eligible": all(gates.values()),
            "source": {
                "embedding": str(candidate_path),
                "embedding_sha256": sha256_path(candidate_path),
                "update_manifest": str(run_path.resolve()),
                "update_manifest_sha256": sha256_path(run_path),
            },
        }

    output = {
        "schema_version": "masld-cl-v13-control-evaluation-v1",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "control_only": True,
        "outcome_variables_read": [],
        "reference_embedding": {
            "path": str(reference_path),
            "sha256": sha256_path(reference_path),
        },
        "harmony_embedding": {
            "path": str(harmony_path),
            "sha256": sha256_path(harmony_path),
        },
        "results": results,
    }
    write_json_exclusive(args.output, output)
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
