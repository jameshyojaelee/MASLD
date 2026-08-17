"""Control-only diagnostic of already trained all-lineage backbones."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _bundle_centroids, _label_scores, _paired_improvement
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .control_evaluation import _neighbor_jaccard_loss, _reference_f1_change_ci
from .embedding import load_embedding, matched_rows
from .firewall import validate_program_firewall


METHODS = {"architecture_surgery", "fine_tune", "replay_only", "ewc_only", "de_novo"}


def _mapping(values: list[str]) -> dict[str, Path]:
    result = {}
    for value in values:
        method, separator, path = value.partition("=")
        if separator != "=" or method in result:
            raise ContractError(f"invalid or duplicate backbone mapping: {value}")
        result[method] = Path(path).resolve()
    if set(result) != METHODS:
        raise ContractError(f"backbone diagnostic requires exactly {sorted(METHODS)}")
    return result


def diagnose_backbones(
    config: dict[str, Any], reference_embedding: str | Path,
    harmony_embedding: str | Path, embedding_values: list[str], output: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    reference_path = Path(reference_embedding).resolve()
    harmony_path = Path(harmony_embedding).resolve()
    reference = load_embedding(reference_path)
    harmony = load_embedding(harmony_path)
    if harmony[0].get("method") != "incumbent_scalesc_pca_harmony":
        raise ContractError("backbone diagnostic Harmony role differs")
    paths = _mapping(embedding_values)
    results = {}
    for index, method in enumerate(sorted(METHODS)):
        path = paths[method]
        bundle = load_embedding(path)
        candidate_reference_positions = np.flatnonzero(
            bundle[2]["strict_reference"].to_numpy(dtype=bool)
        )
        candidate_reference = bundle[2].iloc[candidate_reference_positions].reset_index(drop=True)
        left, right = matched_rows(reference[2], candidate_reference)
        if len(left) != len(reference[2]):
            raise ContractError(f"backbone lacks strict-reference cells: {method}")
        reference_cells = reference[2].iloc[left].reset_index(drop=True)
        candidate_cells = candidate_reference.iloc[right].reset_index(drop=True)
        f1 = _reference_f1_change_ci(
            reference_cells, candidate_cells,
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 1401 + index,
        )
        jaccard = _neighbor_jaccard_loss(
            np.asarray(reference[1])[left],
            np.asarray(bundle[1])[candidate_reference_positions[right]],
            reference_cells,
            k=config["evaluation"]["reference_neighborhood_k"],
            maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
            seed=config["screen"]["seed"] + 1401 + index,
        )
        alignment = _paired_improvement(
            _bundle_centroids(bundle), _bundle_centroids(harmony),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 1501 + index,
        )
        query_macro_f1, query_lineage_f1 = _label_scores(bundle, config["lineages"])
        results[method] = {
            "reference_macro_f1_change": f1,
            "reference_neighborhood_jaccard_loss": jaccard,
            "alignment_vs_harmony": alignment,
            "query_mean_donor_macro_f1": query_macro_f1,
            "query_mean_donor_lineage_f1": query_lineage_f1,
            "gates": {
                "reference_macro_f1": f1["ci_low"] > -config["gates"]["reference_macro_f1_margin"],
                "reference_neighborhood": jaccard <= config["gates"]["reference_neighborhood_jaccard_loss"],
                "alignment": (
                    alignment["improvement"] >= config["gates"]["minimum_control_alignment_improvement"]
                    and alignment["ci_low"] > 0
                ),
            },
            "source": {"path": str(path), "sha256": sha256_path(path)},
        }
    result = {
        "schema_version": "masld-cl-backbone-control-diagnostic-v13",
        "config_sha256": config["_config_sha256"],
        "control_only": True,
        "outcome_variables_read": [],
        "results": results,
        "reference_embedding": {"path": str(reference_path), "sha256": sha256_path(reference_path)},
        "harmony_embedding": {"path": str(harmony_path), "sha256": sha256_path(harmony_path)},
        "use": "prospective V13 architecture diagnosis only; not a promotion result",
    }
    write_json_exclusive(output, result)
    return result
