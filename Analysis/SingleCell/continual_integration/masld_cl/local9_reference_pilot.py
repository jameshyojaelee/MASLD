"""Control-only all-lineage assessment of the nine-donor local reference."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _bundle_centroids, _paired_improvement
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .control_evaluation import _neighbor_jaccard_loss, _reference_f1_change_ci
from .embedding import load_embedding, matched_rows
from .execution import require_execution_ownership, verify_execution_record
from .firewall import validate_program_firewall
from .metrics import macro_f1
from .orthogonal_bridge import (
    _balanced_fit_positions, apply_dataset_control_offsets,
    apply_scaled_orthogonal_bridge, fit_scaled_orthogonal_bridge,
)


def _load_policy(config: dict[str, Any], value: str | Path) -> tuple[Path, dict[str, Any]]:
    pipeline_root = Path(config["_config_path"]).resolve().parent
    path = Path(value).resolve()
    expected = pipeline_root / "reference" / "local9_all_lineage_assessment_policy_v12.json"
    if path != expected:
        raise ContractError("local9 assessment requires the source-controlled V12 policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != "masld-cl-local9-all-lineage-assessment-v12"
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("control_only") is not True
        or policy.get("all_lineage_pilot_only") is not True
        or policy.get("case_stage_program_hero_gene_umap_and_cas13_not_used") is not True
    ):
        raise ContractError("local9 assessment policy is invalid")
    for key in ("training_policy", "reference_assessment"):
        source = pipeline_root / policy[f"{key}_relative_path"]
        if sha256_path(source) != policy[f"{key}_sha256"]:
            raise ContractError(f"local9 assessment policy source changed: {key}")
    return path, policy


def _reclassify(cells, added_donors: set[str]):
    cells = cells.copy()
    added = cells["donor_id"].astype(str).isin(added_donors)
    cells.loc[added, "strict_reference"] = True
    cells.loc[added, "primary_query"] = False
    cells.loc[added, "query_control"] = False
    return cells


def _donor_f1(cells) -> dict[str, float]:
    return {
        str(donor): macro_f1(group["audit_cell_type"], group["predicted_cell_type"])
        for donor, group in cells.groupby("donor_id", sort=True)
    }


def _lineage_counts(cells, donors: set[str]) -> dict[str, dict[str, int]]:
    selected = cells.loc[cells["donor_id"].astype(str).isin(donors)]
    result = {}
    for lineage, group in selected.groupby("audit_cell_type", sort=True):
        result[str(lineage)] = {
            "cells": int(len(group)),
            "donors": int(group["donor_id"].nunique()),
        }
    return result


def evaluate_local9_reference_pilot(
    config: dict[str, Any], policy_value: str | Path,
    current_reference_embedding: str | Path, current_reference_record: str | Path,
    local9_reference_embedding: str | Path, local9_reference_record: str | Path,
    raw_pca_embedding: str | Path, harmony_embedding: str | Path,
    output: str | Path,
) -> dict[str, Any]:
    policy_path, policy = _load_policy(config, policy_value)
    program_firewall_sha256 = validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent

    current_path = Path(current_reference_embedding).resolve()
    local9_path = Path(local9_reference_embedding).resolve()
    current_record_path = Path(current_reference_record).resolve()
    local9_record_path = Path(local9_reference_record).resolve()
    current_record = verify_execution_record(
        current_record_path, pipeline_root, config["_config_sha256"], allow_historical_source=True
    )
    local9_record = verify_execution_record(
        local9_record_path, pipeline_root, config["_config_sha256"], allow_historical_source=True
    )
    require_execution_ownership(current_record, [current_path], role="current seven-donor reference")
    require_execution_ownership(local9_record, [local9_path], role="local9 reference")

    current = load_embedding(current_path)
    local9 = load_embedding(local9_path)
    raw = load_embedding(raw_pca_embedding)
    harmony = load_embedding(harmony_embedding)
    if (
        current[0].get("model_kind") != "all_lineage"
        or local9[0].get("model_kind") != "all_lineage"
        or raw[0].get("method") != "uncorrected_scalesc_pca_30d"
        or harmony[0].get("method") != "incumbent_scalesc_pca_harmony"
    ):
        raise ContractError("local9 assessment received an incompatible embedding role")
    for bundle in (current, local9, raw, harmony):
        if bundle[0].get("config_sha256") != config["_config_sha256"]:
            raise ContractError("local9 assessment embedding config differs")

    added_donors = set(map(str, policy["added_donors"]))
    observed_local9 = set(local9[2]["donor_id"].astype(str))
    if (
        len(observed_local9) != int(policy["challenger_reference_donors"])
        or len(local9[2]) != int(policy["challenger_reference_cells"])
        or not added_donors.issubset(observed_local9)
    ):
        raise ContractError("local9 reference roster differs from policy")

    current_rows, local9_common_rows = matched_rows(current[2], local9[2])
    if len(current_rows) != len(current[2]):
        raise ContractError("local9 reference does not contain the exact current-seven cell roster")
    current_cells = current[2].iloc[current_rows].reset_index(drop=True)
    local9_common_cells = local9[2].iloc[local9_common_rows].reset_index(drop=True)
    f1_change = _reference_f1_change_ci(
        current_cells, local9_common_cells,
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 1201,
    )
    jaccard_loss = _neighbor_jaccard_loss(
        np.asarray(current[1])[current_rows], np.asarray(local9[1])[local9_common_rows], current_cells,
        k=config["evaluation"]["reference_neighborhood_k"],
        maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
        seed=config["screen"]["seed"] + 1201,
    )
    reference_pass = bool(
        f1_change["ci_low"] > -config["gates"]["reference_macro_f1_margin"]
        and jaccard_loss <= config["gates"]["reference_neighborhood_jaccard_loss"]
    )

    raw_rows, local9_rows = matched_rows(raw[2], local9[2])
    order = np.argsort(local9_rows)
    raw_rows, local9_rows = raw_rows[order], local9_rows[order]
    if not np.array_equal(local9_rows, np.arange(len(local9[2]))):
        raise ContractError("raw PCA does not contain the exact ordered local9 roster")
    fit_positions = _balanced_fit_positions(
        local9[2], np.arange(len(local9[2]), dtype=np.int64), "all_lineage", config
    )
    source_mean, target_mean, rotation, scale = fit_scaled_orthogonal_bridge(
        np.asarray(raw[1])[raw_rows[fit_positions]], np.asarray(local9[1])[fit_positions]
    )
    bridged = apply_scaled_orthogonal_bridge(
        np.asarray(raw[1]), source_mean, target_mean, rotation, scale
    )
    bridged[raw_rows] = np.asarray(local9[1], dtype=np.float32)
    candidate_cells = _reclassify(raw[2], added_donors)
    prediction_column = candidate_cells.columns.get_loc("predicted_cell_type")
    candidate_cells.iloc[raw_rows, prediction_column] = (
        local9[2]["predicted_cell_type"].astype(str).to_numpy()
    )

    harmony_cells = _reclassify(harmony[2], added_donors)
    if not np.array_equal(
        candidate_cells["cell_id"].astype(str).to_numpy(),
        harmony_cells["cell_id"].astype(str).to_numpy(),
    ):
        raise ContractError("Harmony and local9 pilot cell orders differ")
    harmony_bundle = (harmony[0], harmony[1], harmony_cells)
    grid = []
    for index, weight in enumerate(map(float, policy["control_offset_weight_grid"])):
        adapted, offsets = apply_dataset_control_offsets(
            bridged, candidate_cells, list(config["roles"]["primary_query_datasets"]),
            weight, int(policy["minimum_control_donors_per_dataset"]),
        )
        comparison = _paired_improvement(
            _bundle_centroids((raw[0], adapted, candidate_cells)),
            _bundle_centroids(harmony_bundle),
            config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 1301 + index,
        )
        alignment_pass = bool(
            comparison["improvement"] >= config["gates"]["minimum_control_alignment_improvement"]
            and comparison["ci_low"] > 0
        )
        grid.append({
            "control_offset_weight": weight,
            "alignment_vs_harmony": comparison,
            "alignment_pass": alignment_pass,
            "reference_pass": reference_pass,
            "eligible": bool(alignment_pass and reference_pass),
            "offsets": offsets,
        })
    eligible = [row for row in grid if row["eligible"]]
    selected = min(
        eligible,
        key=lambda row: (
            row["alignment_vs_harmony"]["candidate_shift"],
            row["control_offset_weight"],
        ),
    ) if eligible else None

    current_f1 = _donor_f1(current[2])
    local9_f1 = _donor_f1(local9[2])
    result = {
        "schema_version": "masld-cl-local9-all-lineage-pilot-v12",
        "config_sha256": config["_config_sha256"],
        "control_only": True,
        "current_reference_is_already_multi_dataset": True,
        "current_reference_summary": {"donors": 7, "datasets": 2, "cells": len(current[2])},
        "challenger_reference_summary": {"donors": 9, "datasets": 3, "cells": len(local9[2])},
        "common_seven_reference": {
            "macro_f1_change": f1_change,
            "neighborhood_jaccard_loss": jaccard_loss,
            "pass": reference_pass,
            "thresholds": {
                "macro_f1_change_ci_low": f"> {-config['gates']['reference_macro_f1_margin']}",
                "neighborhood_jaccard_loss": f"<= {config['gates']['reference_neighborhood_jaccard_loss']}",
            },
        },
        "donor_macro_f1": {"current7": current_f1, "local9": local9_f1},
        "lineage_coverage": {
            "current7": _lineage_counts(local9[2], observed_local9 - added_donors),
            "added2": _lineage_counts(local9[2], added_donors),
        },
        "bridge": {
            "fit_cells": int(len(fit_positions)), "scale": scale,
            "orthonormality_max_abs_error": float(np.max(np.abs(
                rotation.T @ rotation - np.eye(rotation.shape[1])
            ))),
            "diagnostic_only_not_continual_learning": True,
        },
        "control_alignment_grid": grid,
        "selected_control_only_setting": selected,
        "pilot_pass": selected is not None,
        "recommendation": (
            "promising_but_requires_full_identical_comparator_refits"
            if selected is not None else "retain_current_seven_donor_reference"
        ),
        "limitations": [
            "all-lineage pilot only; no lineage-specific challenger models were trained",
            "the closed-form raw-PCA bridge is a diagnostic mapper, not replay-plus-EWC",
            "architecture-surgery and de novo comparators were not refit against local9",
            "no case, stage, program, hero-gene, UMAP, or Cas13 result entered selection",
        ],
        "external_reference_status": {
            "HLiCA_downloaded": True,
            "new_unique_donors": 19,
            "eligible_cells": 142036,
            "mapped_current_4000_hvgs": 3802,
            "missing_current_4000_hvgs": 198,
            "next_step": "requires a prospectively locked common-universe 4000-HVG comparison",
        },
        "program_firewall_sha256": program_firewall_sha256,
        "sources": {
            "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
            "current_reference_embedding": {"path": str(current_path), "sha256": sha256_path(current_path)},
            "current_reference_record": {"path": str(current_record_path), "sha256": sha256_path(current_record_path)},
            "local9_reference_embedding": {"path": str(local9_path), "sha256": sha256_path(local9_path)},
            "local9_reference_record": {"path": str(local9_record_path), "sha256": sha256_path(local9_record_path)},
            "raw_pca_embedding": {"path": str(Path(raw_pca_embedding).resolve()), "sha256": sha256_path(raw_pca_embedding)},
            "harmony_embedding": {"path": str(Path(harmony_embedding).resolve()), "sha256": sha256_path(harmony_embedding)},
        },
    }
    write_json_exclusive(output, result)
    return result
