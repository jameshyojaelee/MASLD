"""Matched-cell disease-preservation evaluation for the locked V9 bridge."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding, matched_rows
from .firewall import load_retargeted_bridge_selection_lock, validate_program_firewall
from .metrics import donor_centroids, donor_distance_spearman, standardized_case_control_separation


def _selected_sources(selection: dict[str, Any]) -> dict[str, dict[str, Any]]:
    selected = float(selection["selected_global_reference_weight"])
    result = {}
    for source in selection["control_results"]:
        with Path(source["path"]).open() as handle:
            row = json.load(handle)
        if float(row["global_reference_weight"]) == selected:
            result[row["model_kind"]] = row
    return result


def _aligned_candidate_and_raw(candidate_manifest: Path, raw_manifest: Path):
    candidate_info, candidate_latent, candidate_cells = load_embedding(candidate_manifest)
    raw_info, raw_latent, raw_cells = load_embedding(raw_manifest)
    if raw_info.get("method") != "uncorrected_scalesc_pca_30d" or raw_info.get("n_latent") != 30:
        raise ContractError("matched outcome comparator is not the frozen uncorrected PC30 authority")
    candidate_rows, raw_rows = matched_rows(candidate_cells, raw_cells)
    if len(candidate_rows) != len(candidate_cells):
        raise ContractError("raw PC30 comparator lacks candidate cells")
    order = np.argsort(candidate_rows)
    candidate_rows, raw_rows = candidate_rows[order], raw_rows[order]
    if not np.array_equal(candidate_rows, np.arange(len(candidate_cells))):
        raise ContractError("candidate/raw matching did not recover the exact candidate order")
    raw_cells = raw_cells.iloc[raw_rows].reset_index(drop=True)
    for column in ("cell_id", "donor_id", "dataset", "audit_cell_type"):
        if not np.array_equal(
            candidate_cells[column].astype(str).to_numpy(),
            raw_cells[column].astype(str).to_numpy(),
        ):
            raise ContractError(f"candidate/raw metadata differs: {column}")
    for column in ("primary_query", "query_control", "analysis_eligible"):
        if not np.array_equal(
            candidate_cells[column].to_numpy(dtype=bool),
            raw_cells[column].to_numpy(dtype=bool),
        ):
            raise ContractError(f"candidate/raw metadata differs: {column}")
    return (
        candidate_info,
        np.asarray(candidate_latent, dtype=np.float64),
        np.asarray(raw_latent[raw_rows], dtype=np.float64),
        candidate_cells,
    )


def _scope_centroids(candidate, raw, cells, lineage: str, study: str | None):
    mask = (
        cells["analysis_eligible"].to_numpy(dtype=bool)
        & cells["primary_query"].to_numpy(dtype=bool)
        & (cells["audit_cell_type"].astype(str).to_numpy() == lineage)
    )
    if study is not None:
        mask &= cells["dataset"].astype(str).to_numpy() == study
    selected = cells.loc[mask].reset_index(drop=True)
    if selected.empty:
        raise ContractError(f"matched outcome scope has no cells: {study}|{lineage}")
    candidate_donors, candidate_centroids = donor_centroids(candidate[mask], selected["donor_id"].astype(str))
    raw_donors, raw_centroids = donor_centroids(raw[mask], selected["donor_id"].astype(str))
    if not np.array_equal(candidate_donors, raw_donors):
        raise ContractError("matched outcome donor rosters differ")
    metadata = selected.drop_duplicates("donor_id").set_index("donor_id").loc[candidate_donors]
    studies = metadata["dataset"].astype(str).to_numpy()
    controls = metadata["query_control"].to_numpy(dtype=bool)
    if study is None:
        for value in sorted(set(studies)):
            positions = studies == value
            candidate_centroids[positions] -= candidate_centroids[positions].mean(axis=0)
            raw_centroids[positions] -= raw_centroids[positions].mean(axis=0)
    return candidate_centroids, raw_centroids, studies, controls, int(mask.sum())


def evaluate_matched_scope(candidate, raw, controls) -> dict[str, float] | None:
    """Compare donor geometry on identical cells; exposed for a numerical fixture."""
    if int(controls.sum()) < 3 or int((~controls).sum()) < 3:
        return None
    candidate_separation = standardized_case_control_separation(candidate[~controls], candidate[controls])
    raw_separation = standardized_case_control_separation(raw[~controls], raw[controls])
    if raw_separation <= 0:
        raise ContractError("matched raw-PC disease separation is zero")
    return {
        "candidate_separation": candidate_separation,
        "raw_pc30_separation": raw_separation,
        "retention": candidate_separation / raw_separation,
        "distance_spearman": donor_distance_spearman(candidate, raw),
    }


def evaluate_orthogonal_bridge_outcomes(
    config: dict[str, Any], selection_lock: str | Path,
    candidate_manifest: str | Path, raw_pca_embedding: str | Path,
    output: str | Path,
) -> list[dict[str, Any]]:
    selection_path = Path(selection_lock).resolve()
    selection = load_retargeted_bridge_selection_lock(selection_path, config)
    validate_program_firewall(config)
    candidate_path = Path(candidate_manifest).resolve()
    raw_path = Path(raw_pca_embedding).resolve()
    candidate_info, candidate, raw, cells = _aligned_candidate_and_raw(candidate_path, raw_path)
    model_kind = candidate_info["model_kind"]
    selected_sources = _selected_sources(selection)
    if set(selected_sources) != {"all_lineage", *config["lineages"]} or model_kind not in selected_sources:
        raise ContractError("locked V9 selection lacks the candidate model")
    selected_result = selected_sources[model_kind]
    if Path(selected_result["sources"]["candidate_embedding"]).resolve() != candidate_path:
        raise ContractError("outcome candidate was not selected by the control-only V9 grid")
    bridge_path = candidate_path.parent / "orthogonal_bridge_manifest.json"
    if sha256_path(bridge_path) != selected_result["sources"]["bridge_manifest_sha256"]:
        raise ContractError("selected V9 bridge manifest changed")
    with bridge_path.open() as handle:
        bridge = json.load(handle)
    parent_path = Path(bridge["parent_manifest_realpath"])
    if sha256_path(parent_path) != bridge["parent_manifest_sha256"]:
        raise ContractError("selected V9 parent bridge changed")
    with parent_path.open() as handle:
        parent = json.load(handle)
    if (
        bridge.get("embedding") != candidate_info
        or float(bridge.get("global_reference_weight", -1)) != float(selection["selected_global_reference_weight"])
        or float(bridge.get("control_offset_weight", -1)) != float(selection["selected_control_offset_weight"])
        or Path(parent.get("raw_pca_embedding_realpath", "")).resolve() != raw_path
        or sha256_path(raw_path) != parent.get("raw_pca_embedding_sha256")
    ):
        raise ContractError("selected V9 outcome provenance is inconsistent")

    lineages = config["lineages"] if model_kind == "all_lineage" else [model_kind]
    rows, scopes = [], []
    for lineage in lineages:
        for study in [None, *config["evaluation"]["powered_query_studies"]]:
            scope = f"POOLED_PRIMARY|{lineage}" if study is None else f"{study}|{lineage}"
            candidate_centroids, raw_centroids, studies, controls, n_cells = _scope_centroids(
                candidate, raw, cells, lineage, study
            )
            result = evaluate_matched_scope(candidate_centroids, raw_centroids, controls)
            detail = {
                "scope": scope, "n_cells": n_cells, "n_donors": int(len(controls)),
                "n_controls": int(controls.sum()), "n_cases": int((~controls).sum()),
                "cell_support": "exactly_matched", "donor_centroid_weighting": "equal_cells_within_donor",
            }
            if result is None:
                scopes.append({**detail, "testable": False, "reason": "fewer_than_three_donors_per_group"})
                continue
            scopes.append({**detail, "testable": True, **result})
            if study is None:
                rows.append({"metric": "pooled_disease_retention", "scope": lineage, "value": result["retention"]})
            else:
                rows.extend([
                    {"metric": "per_study_disease_retention", "scope": scope, "value": result["retention"]},
                    {"metric": "within_study_distance_spearman", "scope": scope, "value": result["distance_spearman"]},
                ])
    if not rows:
        raise ContractError("no V9 matched outcome scope was testable")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json_exclusive(output.with_suffix(output.suffix + ".details.json"), {
        "schema_version": "masld-cl-orthogonal-bridge-matched-outcome-v9",
        "config_sha256": config["_config_sha256"],
        "selection_lock_realpath": str(selection_path),
        "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"],
        "metrics_realpath": str(output.resolve()), "metrics_sha256": sha256_path(output),
        "candidate_embedding": str(candidate_path), "candidate_embedding_sha256": sha256_path(candidate_path),
        "bridge_manifest": str(bridge_path.resolve()), "bridge_manifest_sha256": sha256_path(bridge_path),
        "raw_pca_embedding": str(raw_path), "raw_pca_embedding_sha256": sha256_path(raw_path),
        "model_kind": model_kind,
        "selected_control_offset_weight": selection["selected_control_offset_weight"],
        "selected_global_reference_weight": selection["selected_global_reference_weight"],
        "comparator_contract": "first_30_uncorrected_ScaleSC_PCs_on_the_exact_same_cells_per_donor",
        "original_all_gene_donor_pseudobulk_PCA": "supporting_sensitivity_not_promotion_gate_due_cell_support_mismatch",
        "scopes": scopes, "language_constraint": "cross-sectional stage-associated remodeling",
    })
    return rows
