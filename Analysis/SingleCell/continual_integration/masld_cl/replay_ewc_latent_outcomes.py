"""Post-lock identity and disease-preservation audit for V27."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .firewall import validate_program_firewall
from .latent_knn_label_audit import (
    _fit_predict_knn,
    _prediction_scores,
    load_latent_knn_label_audit_policy,
)
from .orthogonal_bridge_outcome import (
    _aligned_candidate_and_raw,
    _scope_centroids,
    evaluate_matched_scope,
)
from .training import _capped_indices


def evaluate_replay_ewc_latent_adapter_outcomes(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path = Path(policy_value).resolve()
    reference_root = Path(config["_config_path"]).resolve().parent / "reference"
    profiles = {
        reference_root / "replay_ewc_latent_adapter_outcome_policy_v27.json": {
            "policy_schema": "masld-cl-replay-ewc-latent-adapter-outcome-policy-v27",
            "selection_schema": "masld-cl-replay-ewc-latent-adapter-selection-v27",
            "selected_setting": "g100_w100",
            "embedding_method": "reference_frozen_replay_ewc_latent_adapter",
            "result_schema": "masld-cl-replay-ewc-latent-adapter-outcomes-v27",
            "label_source": "final_embedding_knn",
        },
        reference_root / "class_routed_adapter_outcome_policy_v28.json": {
            "policy_schema": "masld-cl-class-routed-adapter-outcome-policy-v28",
            "selection_schema": "masld-cl-class-routed-replay-ewc-adapter-selection-v28",
            "selected_setting": "w075",
            "embedding_method": "class_routed_replay_ewc_geometry_adapter",
            "result_schema": "masld-cl-class-routed-adapter-outcomes-v28",
            "label_source": "final_embedding_knn",
        },
        reference_root / "compartment_routed_adapter_outcome_policy_v29.json": {
            "policy_schema": "masld-cl-compartment-routed-adapter-outcome-policy-v29",
            "selection_schema": "masld-cl-compartment-routed-replay-ewc-adapter-selection-v29",
            "selected_setting": "w075",
            "embedding_method": "compartment_routed_replay_ewc_geometry_adapter",
            "result_schema": "masld-cl-compartment-routed-adapter-outcomes-v29",
            "label_source": "final_embedding_knn",
        },
        reference_root / "immune_compartment_adapter_outcome_policy_v30.json": {
            "policy_schema": "masld-cl-immune-compartment-adapter-outcome-policy-v30",
            "selection_schema": "masld-cl-immune-compartment-replay-ewc-adapter-selection-v30",
            "selected_setting": "w075",
            "embedding_method": "immune_compartment_replay_ewc_geometry_adapter",
            "result_schema": "masld-cl-immune-compartment-adapter-outcomes-v30",
            "label_source": "stored_reference_only_routing_label",
        },
        reference_root / "routed_raw_geometry_outcome_policy_v31.json": {
            "policy_schema": "masld-cl-routed-raw-geometry-outcome-policy-v31",
            "selection_schema": "masld-cl-geometry-preserving-continual-adapter-selection-v14",
            "selected_setting": "g100_w100",
            "embedding_method": "replay_ewc_routed_raw_geometry_adapter",
            "result_schema": "masld-cl-routed-raw-geometry-outcomes-v31",
            "label_source": "stored_reference_only_routing_label",
        },
    }
    profile = profiles.get(policy_path)
    if profile is None:
        raise ContractError("adapter outcomes require a source-controlled policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != profile["policy_schema"]
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("firewall", {}).get("selection_was_locked_before_this_policy") is not True
    ):
        raise ContractError("adapter outcome policy identity differs")
    root = policy_path.parent
    selection_path = (root / policy["selection_lock"]).resolve()
    candidate_path = (root / policy["selected_embedding"]).resolve()
    label_authority_path = (root / policy["label_audit"]["authority"]).resolve()
    raw_path = (root / policy["disease_preservation"]["raw_pca_embedding"]).resolve()
    for path, expected_hash in (
        (selection_path, policy["selection_lock_sha256"]),
        (candidate_path, policy["selected_embedding_sha256"]),
        (label_authority_path, policy["label_audit"]["authority_sha256"]),
        (raw_path, policy["disease_preservation"]["raw_pca_embedding_sha256"]),
    ):
        if sha256_path(path) != expected_hash:
            raise ContractError(f"adapter outcome source changed: {path.name}")
    with selection_path.open() as handle:
        selection = json.load(handle)
    payload = {key: value for key, value in selection.items() if key != "lock_sha256"}
    if (
        selection.get("schema_version") != profile["selection_schema"]
        or selection.get("selection_frozen") is not True
        or selection.get("outcomes_unlocked") is not True
        or selection.get("selected", {}).get("setting") != profile["selected_setting"]
        or selection.get("selected", {}).get("embedding_sha256")
        != policy.get("selection_embedding_sha256", policy["selected_embedding_sha256"])
        or hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
        != selection.get("lock_sha256")
    ):
        raise ContractError("adapter selection lock differs")

    with label_authority_path.open() as handle:
        label_authority = json.load(handle)
    best = policy["label_audit"]["best_non_cl_method_must_be"]
    if (
        label_authority.get("schema_version") != "masld-cl-latent-knn-label-audit-v25"
        or label_authority.get("best_non_cl_method") != best
        or label_authority.get("reference_only_fit") is not True
        or label_authority.get("query_labels_available_to_fit") is not False
    ):
        raise ContractError("adapter label authority differs")
    v25_policy_path, v25_policy = load_latent_knn_label_audit_policy(
        config, label_authority["policy"]["path"]
    )
    if sha256_path(v25_policy_path) != label_authority["policy"]["sha256"]:
        raise ContractError("adapter V25 classifier policy changed")
    candidate_info, candidate_latent, candidate_cells = load_embedding(candidate_path)
    if candidate_info.get("method") != profile["embedding_method"]:
        raise ContractError("selected adapter embedding method differs")
    reference = candidate_cells["strict_reference"].to_numpy(dtype=bool)
    query = candidate_cells["analysis_eligible"].to_numpy(dtype=bool) & ~reference
    if query.sum() != 687559:
        raise ContractError("adapter label-audit roster differs")
    if profile["label_source"] == "stored_reference_only_routing_label":
        if "routing_label" not in candidate_cells:
            raise ContractError("adapter routing-label audit column is missing")
        predictions = candidate_cells.loc[query, "routing_label"].astype(str).to_numpy()
        knn_identity = {
            "source": "stored_reference_only_routing_label",
            "predictions_sha256": hashlib.sha256(
                "\0".join(predictions).encode("utf-8")
            ).hexdigest(),
        }
    else:
        reference_pool = _capped_indices(
            type("ADataView", (), {"obs": candidate_cells, "n_obs": len(candidate_cells)})(),
            np.flatnonzero(reference), "all_lineage", config,
            v25_policy["classifier"]["reference_cap_seed"],
        )
        if len(reference_pool) != 61960:
            raise ContractError("adapter capped reference roster differs")
        predictions, knn_identity = _fit_predict_knn(
            np.asarray(candidate_latent)[reference_pool],
            candidate_cells.iloc[reference_pool]["audit_cell_type"].astype(str).to_numpy(),
            np.asarray(candidate_latent)[query], v25_policy["classifier"],
        )
    candidate_scores = _prediction_scores(
        candidate_cells.loc[query].reset_index(drop=True), predictions, config["lineages"]
    )
    baseline = label_authority["method_scores"][best]
    macro_change = (
        candidate_scores["mean_donor_macro_f1"] - baseline["mean_donor_macro_f1"]
    )
    lineage_changes = {
        lineage: candidate_scores["mean_donor_lineage_f1"][lineage]
        - baseline["mean_donor_lineage_f1"][lineage]
        for lineage in config["lineages"]
    }
    label_gates = {
        "query_macro_f1_noninferior": macro_change
        >= policy["label_audit"]["macro_change_at_least"],
        "all_major_lineages_noninferior": min(lineage_changes.values())
        >= policy["label_audit"]["every_major_lineage_change_at_least"],
    }

    _, candidate, raw, cells = _aligned_candidate_and_raw(candidate_path, raw_path)
    scopes = []
    disease_gates = []
    disease_policy = policy["disease_preservation"]
    for lineage in config["lineages"]:
        for study in [None, *config["evaluation"]["powered_query_studies"]]:
            scope = f"POOLED_PRIMARY|{lineage}" if study is None else f"{study}|{lineage}"
            candidate_centroids, raw_centroids, _, controls, n_cells = _scope_centroids(
                candidate, raw, cells, lineage, study
            )
            matched = evaluate_matched_scope(candidate_centroids, raw_centroids, controls)
            detail = {
                "scope": scope, "n_cells": n_cells, "n_donors": int(len(controls)),
                "n_controls": int(controls.sum()), "n_cases": int((~controls).sum()),
            }
            if matched is None:
                scopes.append({**detail, "testable": False, "reason": "fewer_than_three_donors_per_group"})
                continue
            scopes.append({**detail, "testable": True, **matched})
            if study is None:
                disease_gates.append({
                    "metric": "pooled_disease_retention", "scope": lineage,
                    "observed": matched["retention"],
                    "threshold": disease_policy["pooled_retention_at_least"],
                    "pass": matched["retention"] >= disease_policy["pooled_retention_at_least"],
                })
            else:
                disease_gates.extend(({
                    "metric": "per_study_disease_retention", "scope": scope,
                    "observed": matched["retention"],
                    "threshold": disease_policy["per_study_retention_at_least"],
                    "pass": matched["retention"] >= disease_policy["per_study_retention_at_least"],
                }, {
                    "metric": "within_study_distance_spearman", "scope": scope,
                    "observed": matched["distance_spearman"],
                    "threshold": disease_policy["within_study_distance_spearman_at_least"],
                    "pass": matched["distance_spearman"]
                    >= disease_policy["within_study_distance_spearman_at_least"],
                }))
    result = {
        "schema_version": profile["result_schema"],
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "selection_lock": {"path": str(selection_path), "sha256": sha256_path(selection_path)},
        "candidate_embedding": {"path": str(candidate_path), "sha256": sha256_path(candidate_path)},
        "label_concordance": {
            "classifier": profile["label_source"],
            "query_labels_available_to_fit": False,
            "query_labels_used_only_after_prediction": True,
            "knn_identity": knn_identity,
            "best_non_cl_method": best,
            "candidate_mean_donor_macro_f1": candidate_scores["mean_donor_macro_f1"],
            "baseline_mean_donor_macro_f1": baseline["mean_donor_macro_f1"],
            "query_macro_f1_change": macro_change,
            "candidate_mean_donor_lineage_f1": candidate_scores["mean_donor_lineage_f1"],
            "major_lineage_f1_changes": lineage_changes,
            "gates": label_gates,
            "pass": all(label_gates.values()),
        },
        "disease_preservation": {
            "raw_pca_embedding": {"path": str(raw_path), "sha256": sha256_path(raw_path)},
            "scopes": scopes,
            "gates": disease_gates,
            "pass": bool(disease_gates) and all(row["pass"] for row in disease_gates),
            "language_constraint": disease_policy["language_constraint"],
        },
    }
    result["pilot_pass"] = bool(
        result["label_concordance"]["pass"] and result["disease_preservation"]["pass"]
    )
    write_json_exclusive(output_value, result)
    return result
