"""Control-only selection metrics for common strict7 versus external26."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _bundle_centroids, _label_scores, _paired_improvement
from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path
from .control_evaluation import _neighbor_jaccard_loss, _reference_f1_change_ci
from .embedding import load_embedding, matched_rows
from .external_reference_common import load_external26_policy
from .firewall import validate_program_firewall
from .metrics import macro_f1, positive_class_f1


def _load_owned_reference(config: dict[str, Any], value: str | Path, roster: str):
    path = Path(value).resolve()
    manifest_path = path.parent / "reference_manifest.json"
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    bundle = load_embedding(path)
    if (
        manifest.get("schema_version") != "masld-cl-external26-common-reference-v21"
        or manifest.get("config_sha256") != config["_config_sha256"]
        or manifest.get("roster_name") != roster
        or manifest.get("embedding") != bundle[0]
        or manifest.get("query_labels_hidden") is not True
        or manifest.get("case_stage_program_hero_gene_umap_cas13_used") is not False
    ):
        raise ContractError(f"external26 reference provenance differs: {roster}")
    return bundle, manifest_path, manifest


def _require_identical_common_universe(
    common_manifest: dict[str, Any], expanded_manifest: dict[str, Any],
) -> str:
    common_lock = common_manifest.get("common_prepared_lock", {}).get("lock_sha256")
    expanded_lock = expanded_manifest.get("common_prepared_lock", {}).get("lock_sha256")
    if not common_lock or common_lock != expanded_lock:
        raise ContractError("reference fits do not share one common-universe lock")
    return str(common_lock)


def _held_external_scores(cells, lineages: list[str]) -> dict[str, Any]:
    external = cells.loc[cells["cell_id"].astype(str).str.startswith("HLiCA|")].copy()
    if external.empty or external["donor_id"].nunique() != 19:
        raise ContractError("held-external label evaluation requires 19 donors")
    donor_macro = {
        str(donor): macro_f1(group["audit_cell_type"], group["predicted_cell_type"])
        for donor, group in external.groupby("donor_id", sort=True)
    }
    lineage_scores = {}
    for lineage in lineages:
        values = {}
        for donor, group in external.groupby("donor_id", sort=True):
            truth = group["audit_cell_type"].astype(str).to_numpy() == lineage
            if truth.any():
                prediction = group["predicted_cell_type"].astype(str).to_numpy() == lineage
                values[str(donor)] = positive_class_f1(truth, prediction)
        if not values:
            raise ContractError(f"held external major lineage is not evaluable: {lineage}")
        lineage_scores[lineage] = {
            "donor_positive_f1": values,
            "mean_donor_positive_f1": float(np.mean(list(values.values()))),
            "n_evaluable_donors": len(values),
        }
    return {
        "n_cells": len(external),
        "n_donors": int(external["donor_id"].nunique()),
        "donor_macro_f1": donor_macro,
        "mean_donor_macro_f1": float(np.mean(list(donor_macro.values()))),
        "major_lineages": lineage_scores,
    }


def _external26_gate_decision(
    reference_f1: dict[str, float],
    jaccard_loss: float,
    query_macro_change: float,
    query_lineage_changes: dict[str, float],
    alignment: dict[str, float],
) -> tuple[list[dict[str, Any]], str]:
    """Apply the locked control-only reference expansion gates."""
    worst_lineage = min(query_lineage_changes.values())
    gates = [
        {
            "gate": "current7_macro_f1_change_ci_low",
            "observed": reference_f1["ci_low"], "threshold": -0.02,
            "pass": reference_f1["ci_low"] > -0.02,
        },
        {
            "gate": "current7_neighborhood_jaccard_loss",
            "observed": jaccard_loss, "threshold": 0.05,
            "pass": jaccard_loss <= 0.05,
        },
        {
            "gate": "canonical_query_macro_f1_change",
            "observed": query_macro_change, "threshold": -0.02,
            "pass": query_macro_change >= -0.02,
        },
        {
            "gate": "worst_major_lineage_f1_change",
            "observed": worst_lineage, "threshold": -0.05,
            "pass": worst_lineage >= -0.05,
        },
        {
            "gate": "control_alignment_improvement",
            "observed": alignment["improvement"], "threshold": 0.0,
            "pass": alignment["improvement"] > 0.0 and alignment["ci_low"] > 0.0,
        },
    ]
    decision = (
        "select_external_clean26"
        if all(item["pass"] for item in gates)
        else "retain_common_strict7"
    )
    return gates, decision


def evaluate_external26_references(
    config: dict[str, Any], policy_value: str | Path,
    common_value: str | Path, expanded_value: str | Path,
    output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy = load_external26_policy(config, policy_value)
    common, common_manifest_path, common_manifest = _load_owned_reference(
        config, common_value, "common_strict7"
    )
    expanded, expanded_manifest_path, expanded_manifest = _load_owned_reference(
        config, expanded_value, "external_clean26"
    )
    common_universe_lock = _require_identical_common_universe(
        common_manifest, expanded_manifest
    )
    left, right = matched_rows(common[2], expanded[2])
    if len(left) != len(common[2]) or len(right) != len(expanded[2]):
        raise ContractError("common7 and external26 embeddings have different cell rosters")
    common_latent = np.asarray(common[1])[left]
    common_cells = common[2].iloc[left].reset_index(drop=True)
    expanded_latent = np.asarray(expanded[1])[right]
    expanded_cells = expanded[2].iloc[right].reset_index(drop=True)
    if not np.array_equal(
        common_cells["audit_cell_type"].astype(str).to_numpy(),
        expanded_cells["audit_cell_type"].astype(str).to_numpy(),
    ):
        raise ContractError("external26 frozen audit labels changed")

    current7 = common_cells["strict_reference"].to_numpy(dtype=bool)
    if current7.sum() != 216957 or common_cells.loc[current7, "donor_id"].nunique() != 7:
        raise ContractError("common-universe current7 evaluation roster differs")
    expanded_current7 = expanded_cells.loc[current7].reset_index(drop=True)
    common_current7 = common_cells.loc[current7].reset_index(drop=True)
    reference_f1 = _reference_f1_change_ci(
        common_current7, expanded_current7,
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"],
    )
    jaccard_loss = _neighbor_jaccard_loss(
        common_latent[current7], expanded_latent[current7], common_current7,
        k=config["evaluation"]["reference_neighborhood_k"],
        maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
        seed=config["screen"]["seed"],
    )

    common_query_macro, common_query_lineages = _label_scores(
        (common[0], common_latent, common_cells), config["lineages"]
    )
    expanded_query_macro, expanded_query_lineages = _label_scores(
        (expanded[0], expanded_latent, expanded_cells), config["lineages"]
    )
    query_macro_change = expanded_query_macro - common_query_macro
    query_lineage_changes = {
        lineage: expanded_query_lineages[lineage] - common_query_lineages[lineage]
        for lineage in config["lineages"]
    }
    common_centroids = _bundle_centroids((common[0], common_latent, common_cells))
    expanded_centroids = _bundle_centroids((expanded[0], expanded_latent, expanded_cells))
    alignment = _paired_improvement(
        expanded_centroids, common_centroids,
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 2100,
    )
    held_external = _held_external_scores(common_cells, config["lineages"])

    gates, decision = _external26_gate_decision(
        reference_f1, jaccard_loss, query_macro_change,
        query_lineage_changes, alignment,
    )
    result = {
        "schema_version": "masld-cl-external26-reference-decision-v21",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "sources": {
            "common_strict7_embedding": {"path": str(Path(common_value).resolve()), "sha256": sha256_path(common_value)},
            "common_strict7_manifest": {"path": str(common_manifest_path), "sha256": sha256_path(common_manifest_path)},
            "external_clean26_embedding": {"path": str(Path(expanded_value).resolve()), "sha256": sha256_path(expanded_value)},
            "external_clean26_manifest": {"path": str(expanded_manifest_path), "sha256": sha256_path(expanded_manifest_path)},
        },
        "control_only": True,
        "case_stage_program_hero_gene_umap_cas13_read": False,
        "common_universe_lock_sha256": common_universe_lock,
        "identical_common_universe": True,
        "reference_retention": {
            "macro_f1_change": reference_f1,
            "neighborhood_jaccard_loss": jaccard_loss,
        },
        "query_label_concordance": {
            "common_strict7_macro_f1": common_query_macro,
            "external_clean26_macro_f1": expanded_query_macro,
            "macro_f1_change": query_macro_change,
            "common_strict7_major_lineages": common_query_lineages,
            "external_clean26_major_lineages": expanded_query_lineages,
            "major_lineage_changes": query_lineage_changes,
        },
        "control_alignment": alignment,
        "held_external_transfer_from_common_strict7": held_external,
        "gates": gates,
        "decision": decision,
        "pilot_pass": decision == "select_external_clean26",
        "sensitivity_only": True,
    }
    write_json_exclusive(output_value, result)
    return result


def lock_external26_reference_decision(
    config: dict[str, Any], decision_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    """Freeze the control-only external-reference decision for downstream use."""
    import hashlib

    path = Path(decision_value).resolve()
    with path.open() as handle:
        decision = json.load(handle)
    if (
        decision.get("schema_version")
        != "masld-cl-external26-reference-decision-v21"
        or decision.get("config_sha256") != config["_config_sha256"]
        or decision.get("control_only") is not True
        or decision.get("case_stage_program_hero_gene_umap_cas13_read") is not False
        or decision.get("decision")
        not in {"select_external_clean26", "retain_common_strict7"}
        or decision.get("pilot_pass")
        != (decision.get("decision") == "select_external_clean26")
        or not decision.get("identical_common_universe")
    ):
        raise ContractError("external26 reference decision cannot be locked")
    selected_key = (
        "external_clean26" if decision["pilot_pass"] else "common_strict7"
    )
    source_prefix = selected_key + "_"
    embedding = decision["sources"][source_prefix + "embedding"]
    manifest = decision["sources"][source_prefix + "manifest"]
    for source in (embedding, manifest):
        if sha256_path(source["path"]) != source["sha256"]:
            raise ContractError("selected external26 reference source changed")
    failed_gates = [item["gate"] for item in decision["gates"] if not item["pass"]]
    lock = {
        "schema_version": "masld-cl-external26-reference-lock-v21",
        "config_sha256": config["_config_sha256"],
        "decision": decision["decision"],
        "selected_roster": selected_key,
        "decision_file": {"path": str(path), "sha256": sha256_path(path)},
        "selected_embedding": embedding,
        "selected_manifest": manifest,
        "common_universe_lock_sha256": decision["common_universe_lock_sha256"],
        "failed_gates": failed_gates,
        "case_stage_program_hero_gene_umap_cas13_used": False,
        "expanded_reference_allowed_downstream": bool(decision["pilot_pass"]),
    }
    lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(lock)).hexdigest()
    write_json_exclusive(output_value, lock)
    return lock
