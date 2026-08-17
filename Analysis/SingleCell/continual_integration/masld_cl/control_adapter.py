"""Donor-balanced control calibration of a frozen continual-mapping latent space."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .embedding import load_embedding, matched_rows
from .execution import require_execution_ownership, source_identity, verify_execution_record
from .firewall import validate_program_firewall


def _donor_balanced_centroid(latent: np.ndarray, cells, mask: np.ndarray) -> tuple[np.ndarray, int]:
    positions = np.flatnonzero(mask)
    selected = cells.iloc[positions]
    donor_values = []
    donor_ids = selected["donor_id"].astype(str).to_numpy()
    for donor in sorted(set(donor_ids)):
        donor_values.append(
            np.asarray(latent[positions[donor_ids == donor]], dtype=np.float64).mean(axis=0)
        )
    if not donor_values:
        raise ContractError("control adapter centroid contains no biological donors")
    return np.mean(donor_values, axis=0), len(donor_values)


def fit_and_apply_control_offsets(
    latent: np.ndarray, cells, *, query_datasets: list[str],
    minimum_control_donors: int = 3, minimum_reference_donors: int = 3,
    reference_target_scope: str = "same_preparation",
    global_reference_weight: float = 0.0,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Fit one donor-balanced offset per query technical batch and apply it uniformly."""
    latent = np.asarray(latent)
    if latent.ndim != 2 or len(latent) != len(cells) or not np.isfinite(latent).all():
        raise ContractError("control adapter received an invalid latent matrix")
    if len(set(query_datasets)) != len(query_datasets) or not query_datasets:
        raise ContractError("control adapter query dataset roster is empty or duplicated")
    if reference_target_scope not in {
        "same_preparation", "all_strict_reference", "convex_blend"
    }:
        raise ContractError("control adapter reference target scope is invalid")
    if not np.isfinite(global_reference_weight) or not 0 <= global_reference_weight <= 1:
        raise ContractError("control adapter global reference weight is invalid")
    observed = set(cells.loc[cells["primary_query"], "dataset"].astype(str))
    unknown = sorted(set(query_datasets) - observed)
    if unknown:
        raise ContractError(f"control adapter query datasets are absent: {unknown}")
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    controls = (
        cells["primary_query"].to_numpy(dtype=bool)
        & cells["query_control"].to_numpy(dtype=bool)
        & cells["analysis_eligible"].to_numpy(dtype=bool)
        & ~reference
    )
    requested = cells["dataset"].astype(str).isin(query_datasets).to_numpy()
    batches = sorted(set(cells.loc[controls & requested, "technical_batch"].astype(str)))
    expected_batches = sorted(set(cells.loc[(~reference) & requested, "technical_batch"].astype(str)))
    if batches != expected_batches:
        raise ContractError("every primary query technical batch requires eligible controls")

    result = latent.copy()
    summaries = []
    for batch in batches:
        batch_mask = cells["technical_batch"].astype(str).to_numpy() == batch
        batch_controls = controls & requested & batch_mask
        preparations = sorted(set(cells.loc[batch_controls, "preparation"].astype(str)))
        datasets = sorted(set(cells.loc[batch_controls, "dataset"].astype(str)))
        if len(preparations) != 1 or len(datasets) != 1:
            raise ContractError("technical batch is not one dataset x preparation category")
        preparation = preparations[0]
        control_centroid, n_control = _donor_balanced_centroid(
            latent, cells, batch_controls
        )
        compatible_target = reference & (
            cells["preparation"].astype(str).to_numpy() == preparation
        )
        compatible_centroid, n_compatible = _donor_balanced_centroid(
            latent, cells, compatible_target
        )
        global_centroid, n_global = _donor_balanced_centroid(latent, cells, reference)
        if reference_target_scope == "same_preparation":
            reference_centroid, n_reference = compatible_centroid, n_compatible
        elif reference_target_scope == "all_strict_reference":
            reference_centroid, n_reference = global_centroid, n_global
        else:
            reference_centroid = (
                (1 - global_reference_weight) * compatible_centroid
                + global_reference_weight * global_centroid
            )
            n_reference = n_global
        if (
            n_control < minimum_control_donors
            or n_reference < minimum_reference_donors
            or n_compatible < minimum_reference_donors
        ):
            raise ContractError(
                f"underpowered adapter group {batch}: controls={n_control}, "
                f"reference={n_reference}, compatible_reference={n_compatible}"
            )
        offset = reference_centroid - control_centroid
        apply_mask = (~reference) & requested & batch_mask
        before = np.asarray(latent[apply_mask], dtype=np.float64)
        result[apply_mask] = (
            before + offset
        ).astype(result.dtype, copy=False)
        after = np.asarray(result[apply_mask], dtype=np.float64)
        centered_error = float(np.max(np.abs(
            (after - after.mean(axis=0)) - (before - before.mean(axis=0))
        ), initial=0.0))
        summaries.append({
            "technical_batch": batch,
            "dataset": datasets[0],
            "preparation": preparation,
            "n_control_donors": n_control,
            "n_reference_donors": n_reference,
            "n_compatible_reference_donors": n_compatible,
            "reference_target_scope": reference_target_scope,
            "global_reference_weight": global_reference_weight,
            "n_applied_cells": int(apply_mask.sum()),
            "offset": list(map(float, offset)),
            "offset_l2": float(np.linalg.norm(offset)),
            "maximum_centered_coordinate_error": centered_error,
        })
    if not np.array_equal(result[reference], latent[reference]):
        raise ContractError("control adapter changed strict-reference coordinates")
    return result, summaries


def apply_control_adapter(
    config: dict[str, Any], policy_value: str | Path,
    base_embedding_value: str | Path, base_execution_record_value: str | Path,
    reference_embedding_value: str | Path, reference_execution_record_value: str | Path,
    output_value: str | Path, global_reference_weight: float,
) -> dict[str, Any]:
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    policy_path = Path(policy_value).resolve()
    expected_policy = pipeline_root / "reference" / "control_adapter_policy_v7.json"
    if policy_path != expected_policy:
        raise ContractError("control adapter requires the source-controlled v7 policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != "masld-cl-control-adapter-policy-v7"
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("case_stage_program_hero_gene_and_cas13_outcomes_locked") is not True
    ):
        raise ContractError("control adapter policy is invalid")
    parents = policy.get("parent_control_decisions", [])
    if len(parents) != 6:
        raise ContractError("control adapter v7 requires six parent control decisions")
    parent_results = []
    for source in parents:
        parent = pipeline_root / source["relative_path"]
        if sha256_path(parent) != source.get("sha256"):
            raise ContractError("control adapter parent decision changed")
        with parent.open() as handle:
            parent_results.append(json.load(handle))
    if any(value.get("outcomes_unlocked") is not False for value in parent_results):
        raise ContractError("control adapter parent decisions unlocked outcomes")
    model_kinds = {value.get("model_kind", "all_lineage") for value in parent_results}
    if model_kinds != {"all_lineage", *config["lineages"]}:
        raise ContractError("control adapter parent decisions have the wrong model roster")
    weight = float(global_reference_weight)
    if weight not in policy.get("global_reference_weight_grid", []):
        raise ContractError("control adapter weight is outside the locked v7 grid")

    base_path = Path(base_embedding_value).resolve()
    info, latent, cells = load_embedding(base_path)
    run_path = base_path.parent / "update_manifest.json"
    with run_path.open() as handle:
        run = json.load(handle)
    model_kind = run.get("model_kind")
    if (
        run.get("method") != "architecture_surgery"
        or model_kind not in {"all_lineage", *config["lineages"]}
        or run.get("production") is not False
        or run.get("embedding") != info
        or set(run.get("query_datasets", [])) != set(policy["primary_query_datasets"])
    ):
        raise ContractError("control adapter base is not the locked all-lineage architecture fit")
    record = verify_execution_record(
        base_execution_record_value, pipeline_root, config["_config_sha256"],
        allow_historical_source=True,
    )
    require_execution_ownership(record, [base_path, run_path], role="control adapter base")

    reference_path = Path(reference_embedding_value).resolve()
    reference_info, reference_latent, reference_cells = load_embedding(reference_path)
    reference_run_path = reference_path.parent / "reference_manifest.json"
    with reference_run_path.open() as handle:
        reference_run = json.load(handle)
    if (
        reference_run.get("schema_version") != "masld-cl-reference-v1"
        or reference_run.get("model_kind") != model_kind
        or reference_run.get("embedding") != reference_info
    ):
        raise ContractError("control adapter reference is not the all-lineage frozen reference")
    reference_record = verify_execution_record(
        reference_execution_record_value, pipeline_root, config["_config_sha256"],
        allow_historical_source=True,
    )
    require_execution_ownership(
        reference_record, [reference_path, reference_run_path], role="control adapter reference"
    )
    base_reference_positions = np.flatnonzero(cells["strict_reference"].to_numpy(dtype=bool))
    base_reference_cells = cells.iloc[base_reference_positions].reset_index(drop=True)
    left, right = matched_rows(reference_cells, base_reference_cells)
    if len(left) != len(reference_cells) or len(right) != len(base_reference_cells):
        raise ContractError("frozen reference and architecture reference rosters differ")
    if reference_latent.shape[1] != latent.shape[1]:
        raise ContractError("frozen reference and architecture latent dimensions differ")
    if not np.array_equal(
        reference_cells.iloc[left]["audit_cell_type"].astype(str).to_numpy(),
        base_reference_cells.iloc[right]["audit_cell_type"].astype(str).to_numpy(),
    ):
        raise ContractError("frozen reference audit labels differ from architecture mapping")
    frozen_positions = base_reference_positions[right]
    latent = np.asarray(latent).copy()
    latent[frozen_positions] = np.asarray(reference_latent)[left]
    cells = cells.copy()
    prediction_column = cells.columns.get_loc("predicted_cell_type")
    cells.iloc[frozen_positions, prediction_column] = (
        reference_cells.iloc[left]["predicted_cell_type"].astype(str).to_numpy()
    )
    adapted, offsets = fit_and_apply_control_offsets(
        latent, cells, query_datasets=policy["primary_query_datasets"],
        minimum_control_donors=policy["minimum_control_donors_per_fit_group"],
        minimum_reference_donors=policy["minimum_reference_donors_per_target"],
        reference_target_scope=policy["reference_target_scope"],
        global_reference_weight=weight,
    )
    output = Path(output_value)
    output.mkdir(parents=True, exist_ok=False)
    latent_path = output / "embedding_latent.npy"
    cells_path = output / info["cells_file"]
    np.save(latent_path, adapted)
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(cells.columns)
        writer.writerows(cells.itertuples(index=False, name=None))
    embedding = dict(info)
    embedding.update({
        "latent_file": latent_path.name,
        "cells_file": cells_path.name,
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
        "method": "control_calibrated_continual_adapter",
    })
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    manifest = {
        "schema_version": "masld-cl-control-adapter-v7",
        "config_sha256": config["_config_sha256"],
        "policy_realpath": str(policy_path),
        "policy_sha256": sha256_path(policy_path),
        "base_embedding_realpath": str(base_path),
        "base_embedding_sha256": sha256_path(base_path),
        "base_run_realpath": str(run_path),
        "base_run_sha256": sha256_path(run_path),
        "base_execution_record_realpath": str(Path(base_execution_record_value).resolve()),
        "base_execution_record_sha256": sha256_path(base_execution_record_value),
        "base_source_identity_sha256": record["source_identity_sha256"],
        "reference_embedding_realpath": str(reference_path),
        "reference_embedding_sha256": sha256_path(reference_path),
        "reference_run_realpath": str(reference_run_path),
        "reference_run_sha256": sha256_path(reference_run_path),
        "reference_execution_record_realpath": str(
            Path(reference_execution_record_value).resolve()
        ),
        "reference_execution_record_sha256": sha256_path(reference_execution_record_value),
        "reference_source_identity_sha256": reference_record["source_identity_sha256"],
        "adapter_source_identity": source_identity(pipeline_root),
        "model_kind": model_kind,
        "global_reference_weight": weight,
        "control_only": True,
        "outcomes_unlocked": False,
        "reference_coordinates_and_predictions_frozen_to_preupdate": True,
        "offsets": offsets,
        "embedding": embedding,
    }
    write_json_exclusive(output / "control_adapter_manifest.json", manifest)
    return manifest
