"""Diagnose the cell-centroid versus donor-pseudobulk preservation comparator."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py

from .config import repo_path, write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .firewall import load_control_adapter_selection_lock, validate_program_firewall
from .outcome_evaluation import _candidate_donors, _evaluate_scope, _load_raw_scope, _match_raw


def diagnose_outcome_comparator(
    config: dict[str, Any], selection_lock: str | Path,
    harmony_manifest: str | Path, raw_pca_manifest: str | Path,
    output: str | Path,
) -> dict[str, Any]:
    selection_path = Path(selection_lock).resolve()
    selection = load_control_adapter_selection_lock(selection_path, config)
    validate_program_firewall(config)
    harmony_path = Path(harmony_manifest).resolve()
    harmony_info, harmony, cells = load_embedding(harmony_path)
    if harmony_info.get("method") != "incumbent_scalesc_pca_harmony":
        raise ContractError("comparator diagnostic requires the audited Harmony bundle")
    atlas = repo_path(config, config["input"]["atlas_h5ad"])
    with h5py.File(atlas, "r") as handle:
        uncorrected = handle["obsm/X_pca"][:]
    if uncorrected.shape[0] != len(cells):
        raise ContractError("uncorrected PCA and Harmony cell rosters differ")
    raw_path = Path(raw_pca_manifest).resolve()
    rows = []
    for method, latent in (("harmony", harmony), ("uncorrected_cell_pca", uncorrected)):
        for lineage in config["lineages"]:
            for study in config["evaluation"]["powered_query_studies"]:
                scope = f"{study}|{lineage}"
                manifest, _, raw = _load_raw_scope(raw_path, scope)
                if manifest.get("selection_lock_sha256") != selection["lock_sha256"]:
                    raise ContractError("raw-PCA authority differs from v7 selection lock")
                donors, candidate, studies, controls = _candidate_donors(
                    cells, latent, lineage, study
                )
                raw_latent = _match_raw(donors, candidate, studies, controls, raw)
                result = _evaluate_scope(candidate, raw_latent, controls)
                if result is None:
                    raise ContractError(f"unexpected untestable powered scope: {scope}")
                rows.append({"method": method, "scope": scope, **result})
    threshold = float(config["gates"]["within_study_distance_spearman"])
    result = {
        "schema_version": "masld-cl-outcome-comparator-diagnostic-v1",
        "config_sha256": config["_config_sha256"],
        "selection_lock": str(selection_path),
        "selection_lock_file_sha256": sha256_path(selection_path),
        "raw_pca_manifest": str(raw_path),
        "raw_pca_manifest_sha256": sha256_path(raw_path),
        "harmony_manifest": str(harmony_path),
        "harmony_manifest_sha256": sha256_path(harmony_path),
        "atlas": str(atlas.resolve()),
        "atlas_sha256": sha256_path(atlas),
        "distance_spearman_threshold": threshold,
        "rows": rows,
        "summary": {
            method: {
                "minimum_distance_spearman": min(
                    row["distance_spearman"] for row in rows if row["method"] == method
                ),
                "maximum_distance_spearman": max(
                    row["distance_spearman"] for row in rows if row["method"] == method
                ),
                "passing_scopes": sum(
                    row["distance_spearman"] >= threshold
                    for row in rows if row["method"] == method
                ),
                "total_scopes": sum(row["method"] == method for row in rows),
            }
            for method in ("harmony", "uncorrected_cell_pca")
        },
        "interpretation": (
            "The prespecified comparison mixes equal-cell latent centroids with "
            "library-size-weighted donor-pseudobulk log-CPM PCA; failure of both the "
            "incumbent and uncorrected cell PCA indicates a comparator mismatch, not "
            "evidence that donor-aware embedding correction is permissible."
        ),
    }
    write_json_exclusive(output, result)
    return result
