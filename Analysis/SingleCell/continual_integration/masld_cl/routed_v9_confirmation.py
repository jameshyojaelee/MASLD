"""Locked V32 label, seed, order, and fixed-weight held-study confirmation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .contracts import ContractError, sha256_path
from .embedding import load_embedding
from .firewall import validate_program_firewall
from .latent_knn_label_audit import _prediction_scores
from .orthogonal_bridge_held_study import _grid, _study_result, _summarize
from .config import write_json_exclusive


def confirm_routed_v9(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path = Path(policy_value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "routed_v9_confirmation_policy_v32.json"
    )
    if policy_path != expected:
        raise ContractError("V32 confirmation requires its source-controlled policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != "masld-cl-routed-v9-confirmation-policy-v32"
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("held_study_rule", {}).get("reselect_weight_inside_training_study") is not False
    ):
        raise ContractError("V32 confirmation policy identity differs")
    root = policy_path.parent
    fields = (
        "candidate_embedding", "selection_lock", "outcome_decision",
        "seed_confirmation", "order_confirmation", "label_authority",
    )
    paths = {}
    for field in fields:
        path = (root / policy[field]).resolve()
        if sha256_path(path) != policy[f"{field}_sha256"]:
            raise ContractError(f"V32 confirmation source changed: {field}")
        paths[field] = path
    with paths["selection_lock"].open() as handle:
        selection = json.load(handle)
    weight = float(policy["global_reference_weight"])
    summary = selection.get("summaries", {}).get(str(weight))
    if (
        float(selection.get("selected_global_reference_weight", -1)) != weight
        or not summary or summary.get("eligible") is not True
        or summary.get("preservation") is not True
    ):
        raise ContractError("V32 source selection is not the locked eligible g075 bundle")
    with paths["outcome_decision"].open() as handle:
        outcomes = json.load(handle)
    with paths["seed_confirmation"].open() as handle:
        seeds = json.load(handle)
    with paths["order_confirmation"].open() as handle:
        orders = json.load(handle)
    if (
        outcomes.get("disease_preservation_pass") is not True
        or any(model.get("pass") is not True for model in outcomes.get("models", {}).values())
        or seeds.get("all_five_seeds_pass") is not True
        or orders.get("order_gate_pass") is not True
        or orders.get("all_six_orders_coordinate_equivalent") is not True
    ):
        raise ContractError("V32 inherited disease, seed, or order authority did not pass")

    _, _, cells = load_embedding(paths["candidate_embedding"])
    query = cells["analysis_eligible"].to_numpy(dtype=bool) & ~cells["strict_reference"].to_numpy(dtype=bool)
    if "routing_label" not in cells or query.sum() != 687559:
        raise ContractError("V32 routing-label roster differs")
    predictions = cells.loc[query, "routing_label"].astype(str).to_numpy()
    candidate_scores = _prediction_scores(
        cells.loc[query].reset_index(drop=True), predictions, config["lineages"]
    )
    with paths["label_authority"].open() as handle:
        label_authority = json.load(handle)
    best = policy["best_non_cl_method_must_be"]
    if label_authority.get("best_non_cl_method") != best:
        raise ContractError("V32 best non-CL label authority differs")
    baseline = label_authority["method_scores"][best]
    macro_change = candidate_scores["mean_donor_macro_f1"] - baseline["mean_donor_macro_f1"]
    lineage_changes = {
        lineage: candidate_scores["mean_donor_lineage_f1"][lineage]
        - baseline["mean_donor_lineage_f1"][lineage]
        for lineage in config["lineages"]
    }
    label_gates = {
        "query_macro_f1_noninferior": macro_change >= policy["gates"]["macro_change_at_least"],
        "all_major_lineages_noninferior": min(lineage_changes.values())
        >= policy["gates"]["every_major_lineage_change_at_least"],
    }

    grid = _grid(selection)
    kinds = ("all_lineage", *config["lineages"])
    cache = {}
    held = {}
    for study_index, study in enumerate(policy["held_study_rule"]["studies"]):
        rows = [
            _study_result(
                config, grid[(weight, kind)], study, cache,
                15001 + 100 * study_index + 10 * kind_index,
            )
            for kind_index, kind in enumerate(kinds)
        ]
        held[study] = {
            "global_reference_weight": weight,
            "reselected_inside_training_study": False,
            "rows": rows,
            "summary": _summarize(config, rows),
            "query_cell_type_and_stage_labels_hidden": True,
            "held_control_indicator_used_only_for_control_calibration": True,
            "disease_geometry_pass_inherited_from_uniform_within_study_translation": True,
        }
    held_pass = all(value["summary"]["eligible"] for value in held.values())
    result = {
        "schema_version": "masld-cl-routed-v9-confirmation-v32",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "control_selection": summary,
        "disease_preservation_pass": True,
        "label_concordance": {
            "source": "stored_reference_only_replay_ewc_routing_label",
            "query_labels_available_to_fit": False,
            "query_macro_f1_change": macro_change,
            "major_lineage_f1_changes": lineage_changes,
            "gates": label_gates,
            "pass": all(label_gates.values()),
        },
        "seed_confirmation": {
            "all_five_seeds_pass": True,
            "source": str(paths["seed_confirmation"]),
        },
        "order_confirmation": {
            "all_six_orders_pass": True,
            "minimum_distance_spearman": orders["minimum_order_distance_spearman"],
            "worst_order_degradation": orders["worst_order_degradation"],
            "order_to_de_novo_seed_variability_ratio": orders[
                "order_to_de_novo_seed_variability_ratio"
            ],
        },
        "held_studies": held,
        "fixed_weight_held_studies_pass": held_pass,
    }
    result["confirmation_pass"] = bool(
        result["label_concordance"]["pass"]
        and result["disease_preservation_pass"]
        and result["seed_confirmation"]["all_five_seeds_pass"]
        and result["order_confirmation"]["all_six_orders_pass"]
        and held_pass
    )
    write_json_exclusive(output_value, result)
    return result
