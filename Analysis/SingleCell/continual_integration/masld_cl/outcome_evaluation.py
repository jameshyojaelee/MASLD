"""Post-selection disease-state preservation metrics at donor level."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .execution import require_execution_ownership, verify_execution_record
from .firewall import (
    load_control_adapter_selection_lock,
    load_selection_lock,
    validate_program_firewall,
)
from .metrics import donor_centroids, donor_distance_spearman, standardized_case_control_separation


def _load_raw_scope(manifest_path: Path, scope: str):
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    if manifest.get("schema_version") != "masld-cl-raw-pca-v1":
        raise ContractError("unsupported raw-PCA manifest")
    entries = {entry["scope"]: entry for entry in manifest["scopes"]}
    if scope not in entries:
        raise ContractError(f"raw-PCA scope is missing: {scope}")
    entry = entries[scope]
    path = manifest_path.parent / entry["file"]
    if sha256_path(path) != entry["sha256"]:
        raise ContractError(f"raw-PCA artifact changed: {scope}")
    data = np.load(path)
    return manifest, entry, data


def _candidate_donors(cells, latent, lineage: str, study: str | None):
    mask = (
        cells["analysis_eligible"].to_numpy(dtype=bool)
        & cells["primary_query"].to_numpy(dtype=bool)
        & (cells["audit_cell_type"].astype(str).to_numpy() == lineage)
    )
    if study is not None:
        mask &= cells["dataset"].astype(str).to_numpy() == study
    selected = cells.loc[mask].reset_index(drop=True)
    donors, centroids = donor_centroids(
        np.asarray(latent[mask], dtype=np.float64), selected["donor_id"].astype(str)
    )
    metadata = selected.drop_duplicates("donor_id").set_index("donor_id").loc[donors]
    studies = metadata["dataset"].astype(str).to_numpy()
    controls = metadata["query_control"].to_numpy(dtype=bool)
    if study is None:
        for value in sorted(set(studies)):
            centroids[studies == value] -= centroids[studies == value].mean(axis=0)
    return donors.astype(str), centroids, studies, controls


def _match_raw(candidate_donors, candidate, studies, controls, raw):
    raw_donors = raw["donor_id"].astype(str)
    raw_index = {donor: i for i, donor in enumerate(raw_donors)}
    if set(candidate_donors) != set(raw_donors):
        missing = sorted(set(raw_donors) - set(candidate_donors))
        extra = sorted(set(candidate_donors) - set(raw_donors))
        raise ContractError(f"candidate/raw donor mismatch; missing={missing}, extra={extra}")
    order = np.asarray([raw_index[x] for x in candidate_donors], dtype=np.int64)
    raw_control = raw["query_control"][order].astype(bool)
    raw_study = raw["study"][order].astype(str)
    if not np.array_equal(raw_control, controls) or not np.array_equal(raw_study, studies):
        raise ContractError("candidate/raw donor metadata mismatch")
    return np.asarray(raw["latent"][order], dtype=np.float64)


def _evaluate_scope(candidate, raw, controls):
    if controls.sum() < 3 or (~controls).sum() < 3:
        return None
    candidate_separation = standardized_case_control_separation(
        candidate[~controls], candidate[controls]
    )
    raw_separation = standardized_case_control_separation(raw[~controls], raw[controls])
    if raw_separation <= 0:
        raise ContractError("raw-PCA disease separation is zero")
    return {
        "candidate_separation": candidate_separation,
        "raw_pca_separation": raw_separation,
        "retention": candidate_separation / raw_separation,
        "distance_spearman": donor_distance_spearman(candidate, raw),
    }


def evaluate_outcome_preservation(
    config: dict[str, Any], selection_lock: str | Path,
    candidate_manifest: str | Path, raw_pca_manifest: str | Path,
    execution_record: str | Path, output: str | Path,
) -> list[dict[str, Any]]:
    selection = load_selection_lock(selection_lock, config)
    validate_program_firewall(config)
    candidate_info, latent, cells = load_embedding(candidate_manifest)
    candidate_manifest = Path(candidate_manifest).resolve()
    candidate_run_path = candidate_manifest.parent / "update_manifest.json"
    if not candidate_run_path.is_file():
        raise ContractError("outcome candidate lacks sibling update manifest")
    with candidate_run_path.open() as handle:
        candidate_run = json.load(handle)
    if (
        candidate_run.get("config_sha256") != config["_config_sha256"]
        or candidate_run.get("selection_lock_sha256") != selection["lock_sha256"]
        or candidate_run.get("embedding") != candidate_info
        or candidate_run.get("method") != "continual_learning"
        or candidate_run.get("replay_mode") != "random"
        or candidate_run.get("sensitivity_only") is True
        or float(candidate_run.get("ewc_lambda")) != float(selection["selected"]["ewc_lambda"])
        or float(candidate_run.get("replay_fraction")) != float(selection["selected"]["replay_fraction"])
    ):
        raise ContractError("outcome candidate provenance differs from locked selection")
    pipeline_root = Path(config["_config_path"]).resolve().parent
    record_path = Path(execution_record).resolve()
    record = verify_execution_record(record_path, pipeline_root, config["_config_sha256"])
    require_execution_ownership(
        record, [candidate_manifest, candidate_run_path], role="outcome candidate"
    )
    if candidate_info["model_kind"] == "all_lineage":
        lineages = config["lineages"]
    else:
        lineages = [candidate_info["model_kind"]]
    rows, details, raw_scope_sources = [], [], []
    raw_pca_manifest = Path(raw_pca_manifest).resolve()
    for lineage in lineages:
        for study in [None, *config["evaluation"]["powered_query_studies"]]:
            scope = f"POOLED_PRIMARY|{lineage}" if study is None else f"{study}|{lineage}"
            manifest, entry, raw = _load_raw_scope(raw_pca_manifest, scope)
            raw_scope_path = raw_pca_manifest.parent / entry["file"]
            raw_scope_sources.append({
                "scope": scope, "path": str(raw_scope_path.resolve()),
                "sha256": entry["sha256"],
            })
            if manifest.get("selection_lock_sha256") != selection["lock_sha256"]:
                raise ContractError("raw-PCA authority predates or differs from selection lock")
            donors, candidate, studies, controls = _candidate_donors(
                cells, latent, lineage, study
            )
            raw_latent = _match_raw(donors, candidate, studies, controls, raw)
            result = _evaluate_scope(candidate, raw_latent, controls)
            if result is None:
                details.append({
                    "scope": scope, "testable": False, "n_controls": int(controls.sum()),
                    "n_cases": int((~controls).sum()), "reason": "fewer_than_three_donors_per_group",
                })
                continue
            details.append({
                "scope": scope, "testable": True, "n_controls": int(controls.sum()),
                "n_cases": int((~controls).sum()), **result,
            })
            if study is None:
                rows.append({"metric": "pooled_disease_retention", "scope": lineage, "value": result["retention"]})
            else:
                rows.extend([
                    {"metric": "per_study_disease_retention", "scope": scope, "value": result["retention"]},
                    {"metric": "within_study_distance_spearman", "scope": scope, "value": result["distance_spearman"]},
                ])
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ContractError("no powered outcome-preservation scope was testable")
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json_exclusive(output.with_suffix(output.suffix + ".details.json"), {
        "schema_version": "masld-cl-outcome-evaluation-v1",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "metrics_realpath": str(output.resolve()),
        "metrics_sha256": sha256_path(output),
        "candidate_embedding": str(candidate_manifest),
        "candidate_embedding_sha256": sha256_path(candidate_manifest),
        "candidate_run_manifest": str(candidate_run_path.resolve()),
        "candidate_run_manifest_sha256": sha256_path(candidate_run_path),
        "execution_records": [{"path": str(record_path), "sha256": sha256_path(record_path)}],
        "raw_pca_manifest": str(raw_pca_manifest),
        "raw_pca_manifest_sha256": sha256_path(raw_pca_manifest),
        "raw_pca_scope_sources": raw_scope_sources,
        "model_kind": candidate_info["model_kind"],
        "scopes": details,
        "language_constraint": "cross-sectional stage-associated remodeling",
    })
    return rows


def evaluate_control_adapter_outcome_preservation(
    config: dict[str, Any], selection_lock: str | Path,
    candidate_manifest: str | Path, raw_pca_manifest: str | Path,
    output: str | Path,
) -> list[dict[str, Any]]:
    """Evaluate v7 after its control-only lock, preserving the original EWC path."""
    selection_path = Path(selection_lock).resolve()
    selection = load_control_adapter_selection_lock(selection_path, config)
    validate_program_firewall(config)
    candidate_info, latent, cells = load_embedding(candidate_manifest)
    candidate_manifest = Path(candidate_manifest).resolve()
    adapter_path = candidate_manifest.parent / "control_adapter_manifest.json"
    if not adapter_path.is_file():
        raise ContractError("control-adapter outcome candidate lacks its adapter manifest")
    with adapter_path.open() as handle:
        adapter = json.load(handle)
    selected_weight = float(selection["selected_global_reference_weight"])
    expected_kind = candidate_info["model_kind"]
    if (
        adapter.get("schema_version") != "masld-cl-control-adapter-v7"
        or adapter.get("config_sha256") != config["_config_sha256"]
        or adapter.get("embedding") != candidate_info
        or adapter.get("model_kind") != expected_kind
        or float(adapter.get("global_reference_weight", -1)) != selected_weight
        or adapter.get("control_only") is not True
        or adapter.get("outcomes_unlocked") is not False
        or adapter.get("adapter_source_identity") != selection.get("source_identity")
    ):
        raise ContractError("control-adapter outcome candidate differs from locked selection")
    selected_manifests = set()
    for source in selection["control_results"]:
        with Path(source["path"]).open() as handle:
            result = json.load(handle)
        if float(result["global_reference_weight"]) == selected_weight:
            selected_manifests.add(str(Path(result["sources"]["adapter_manifest"]).resolve()))
    if str(adapter_path) not in selected_manifests:
        raise ContractError("outcome candidate was not selected by the control-only grid")
    for path_key, hash_key in (
        ("base_embedding_realpath", "base_embedding_sha256"),
        ("base_execution_record_realpath", "base_execution_record_sha256"),
        ("base_run_realpath", "base_run_sha256"),
        ("reference_embedding_realpath", "reference_embedding_sha256"),
        ("reference_execution_record_realpath", "reference_execution_record_sha256"),
        ("reference_run_realpath", "reference_run_sha256"),
        ("policy_realpath", "policy_sha256"),
    ):
        if sha256_path(adapter.get(path_key, "")) != adapter.get(hash_key):
            raise ContractError(f"control-adapter provenance changed: {path_key}")
    if candidate_info["model_kind"] == "all_lineage":
        lineages = config["lineages"]
    else:
        lineages = [candidate_info["model_kind"]]
    rows, details, raw_scope_sources = [], [], []
    raw_pca_manifest = Path(raw_pca_manifest).resolve()
    for lineage in lineages:
        for study in [None, *config["evaluation"]["powered_query_studies"]]:
            scope = f"POOLED_PRIMARY|{lineage}" if study is None else f"{study}|{lineage}"
            manifest, entry, raw = _load_raw_scope(raw_pca_manifest, scope)
            raw_scope_path = raw_pca_manifest.parent / entry["file"]
            raw_scope_sources.append({
                "scope": scope, "path": str(raw_scope_path.resolve()),
                "sha256": entry["sha256"],
            })
            if manifest.get("selection_lock_sha256") != selection["lock_sha256"]:
                raise ContractError("raw-PCA authority differs from v7 selection lock")
            donors, candidate, studies, controls = _candidate_donors(
                cells, latent, lineage, study
            )
            raw_latent = _match_raw(donors, candidate, studies, controls, raw)
            result = _evaluate_scope(candidate, raw_latent, controls)
            if result is None:
                details.append({
                    "scope": scope, "testable": False, "n_controls": int(controls.sum()),
                    "n_cases": int((~controls).sum()), "reason": "fewer_than_three_donors_per_group",
                })
                continue
            details.append({
                "scope": scope, "testable": True, "n_controls": int(controls.sum()),
                "n_cases": int((~controls).sum()), **result,
            })
            if study is None:
                rows.append({"metric": "pooled_disease_retention", "scope": lineage, "value": result["retention"]})
            else:
                rows.extend([
                    {"metric": "per_study_disease_retention", "scope": scope, "value": result["retention"]},
                    {"metric": "within_study_distance_spearman", "scope": scope, "value": result["distance_spearman"]},
                ])
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ContractError("no powered control-adapter outcome scope was testable")
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "scope", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json_exclusive(output.with_suffix(output.suffix + ".details.json"), {
        "schema_version": "masld-cl-control-adapter-outcome-v7",
        "config_sha256": config["_config_sha256"],
        "selection_lock_realpath": str(selection_path),
        "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"],
        "metrics_realpath": str(output.resolve()),
        "metrics_sha256": sha256_path(output),
        "candidate_embedding": str(candidate_manifest),
        "candidate_embedding_sha256": sha256_path(candidate_manifest),
        "adapter_manifest": str(adapter_path),
        "adapter_manifest_sha256": sha256_path(adapter_path),
        "raw_pca_manifest": str(raw_pca_manifest),
        "raw_pca_manifest_sha256": sha256_path(raw_pca_manifest),
        "raw_pca_scope_sources": raw_scope_sources,
        "model_kind": candidate_info["model_kind"],
        "global_reference_weight": selected_weight,
        "scopes": details,
        "language_constraint": "cross-sectional stage-associated remodeling",
    })
    return rows
