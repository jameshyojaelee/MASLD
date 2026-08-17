"""Control-only all-lineage pilot gate before the full six-model screen."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import (
    _bundle_centroids, _paired_improvement, _validate_harmony_info,
)
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .control_evaluation import _neighbor_jaccard_loss, _reference_f1_change_ci
from .embedding import load_embedding, matched_rows
from .execution import require_execution_ownership, verify_execution_record
from .firewall import validate_program_firewall


def _pilot_passes(config: dict[str, Any], metrics: dict[str, float]) -> dict[str, bool]:
    gates = config["gates"]
    return {
        "reference_retention": (
            metrics["reference_macro_f1_change_ci_low"]
            > -gates["reference_macro_f1_margin"]
            and metrics["reference_neighborhood_jaccard_loss"]
            <= gates["reference_neighborhood_jaccard_loss"]
        ),
        "no_catastrophic_control_worsening": (
            metrics["alignment_improvement_vs_harmony"]
            >= -gates["maximum_stratum_worsening"]
            and metrics["alignment_improvement_vs_architecture_surgery"]
            >= -gates["maximum_stratum_worsening"]
        ),
    }


def _run_identity(
    config: dict[str, Any], embedding_value: str | Path,
    record_value: str | Path, expected_schema: str,
) -> tuple[tuple[dict[str, Any], np.ndarray, Any], dict[str, Any], dict[str, Any], Path]:
    embedding_path = Path(embedding_value).resolve()
    bundle = load_embedding(embedding_path)
    manifest_name = (
        "reference_manifest.json" if expected_schema == "masld-cl-reference-v1"
        else "update_manifest.json"
    )
    run_path = embedding_path.parent / manifest_name
    with run_path.open() as handle:
        run = json.load(handle)
    if (
        run.get("schema_version") != expected_schema
        or run.get("config_sha256") != config["_config_sha256"]
        or run.get("model_kind") != "all_lineage"
        or run.get("embedding") != bundle[0]
    ):
        raise ContractError("pilot training manifest does not own its all-lineage embedding")
    record_path = Path(record_value).resolve()
    record = verify_execution_record(
        record_path, Path(config["_config_path"]).resolve().parent,
        config["_config_sha256"],
    )
    require_execution_ownership(record, [embedding_path, run_path], role="pilot")
    return bundle, run, record, run_path


def evaluate_all_lineage_pilot(
    config: dict[str, Any], reference_embedding: str | Path,
    harmony_embedding: str | Path, candidate_embedding: str | Path,
    architecture_embedding: str | Path, reference_execution_record: str | Path,
    candidate_execution_record: str | Path,
    architecture_execution_record: str | Path, output: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    reference, reference_run, _, reference_run_path = _run_identity(
        config, reference_embedding, reference_execution_record,
        "masld-cl-reference-v1",
    )
    candidate, candidate_run, _, candidate_run_path = _run_identity(
        config, candidate_embedding, candidate_execution_record,
        "masld-cl-update-v1",
    )
    architecture, architecture_run, _, architecture_run_path = _run_identity(
        config, architecture_embedding, architecture_execution_record,
        "masld-cl-update-v1",
    )
    expected_query = set(config["evaluation"]["powered_query_studies"])
    if (
        candidate_run.get("method") != "continual_learning"
        or candidate_run.get("replay_mode") != "random"
        or candidate_run.get("sensitivity_only") is not False
        or candidate_run.get("production") is not False
        or candidate_run.get("seed") != config["screen"]["seed"]
        or float(candidate_run.get("ewc_lambda")) != 100.0
        or float(candidate_run.get("replay_fraction")) != 0.20
        or set(candidate_run.get("query_datasets", [])) != expected_query
        or candidate_run.get("held_out_datasets")
    ):
        raise ContractError("pilot candidate is not the seed-17 paper setting")
    if (
        architecture_run.get("method") != "architecture_surgery"
        or architecture_run.get("production") is not False
        or architecture_run.get("seed") != config["screen"]["seed"]
        or float(architecture_run.get("ewc_lambda")) != 0
        or float(architecture_run.get("replay_fraction")) != 0
        or set(architecture_run.get("query_datasets", [])) != expected_query
        or architecture_run.get("held_out_datasets")
    ):
        raise ContractError("pilot architecture comparator has the wrong context")
    harmony = load_embedding(harmony_embedding)
    _validate_harmony_info(harmony[0], config)
    reference_cells = reference[2].reset_index(drop=True)
    candidate_reference = candidate[2].loc[
        candidate[2]["strict_reference"]
    ].reset_index(drop=True)
    left, right = matched_rows(reference_cells, candidate_reference)
    f1 = _reference_f1_change_ci(
        reference_cells.iloc[left].reset_index(drop=True),
        candidate_reference.iloc[right].reset_index(drop=True),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"],
    )
    candidate_reference_indices = np.flatnonzero(
        candidate[2]["strict_reference"].to_numpy(dtype=bool)
    )[right]
    jaccard = _neighbor_jaccard_loss(
        np.asarray(reference[1])[left],
        np.asarray(candidate[1])[candidate_reference_indices],
        reference_cells.iloc[left].reset_index(drop=True),
        k=config["evaluation"]["reference_neighborhood_k"],
        maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
        seed=config["screen"]["seed"],
    )
    candidate_centroids = _bundle_centroids(candidate)
    harmony_result = _paired_improvement(
        candidate_centroids, _bundle_centroids(harmony),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 701,
    )
    architecture_result = _paired_improvement(
        candidate_centroids, _bundle_centroids(architecture),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 801,
    )
    metrics = {
        "reference_macro_f1_change_ci_low": f1["ci_low"],
        "reference_neighborhood_jaccard_loss": jaccard,
        "alignment_improvement_vs_harmony": harmony_result["improvement"],
        "alignment_improvement_vs_architecture_surgery": architecture_result["improvement"],
    }
    gates = _pilot_passes(config, metrics)
    result = {
        "schema_version": "masld-cl-all-lineage-pilot-v1",
        "config_sha256": config["_config_sha256"],
        "metrics": metrics,
        "bootstrap": {
            "harmony": harmony_result,
            "architecture_surgery": architecture_result,
        },
        "gates": gates,
        "passed": all(gates.values()),
        "control_only": True,
        "sources": {
            "reference_run": {"path": str(reference_run_path), "sha256": sha256_path(reference_run_path)},
            "candidate_run": {"path": str(candidate_run_path), "sha256": sha256_path(candidate_run_path)},
            "architecture_run": {"path": str(architecture_run_path), "sha256": sha256_path(architecture_run_path)},
            "harmony_embedding": {"path": str(Path(harmony_embedding).resolve()), "sha256": sha256_path(harmony_embedding)},
            "execution_records": [
                {"path": str(Path(value).resolve()), "sha256": sha256_path(value)}
                for value in (
                    reference_execution_record, candidate_execution_record,
                    architecture_execution_record,
                )
            ],
        },
    }
    write_json_exclusive(output, result)
    return result
