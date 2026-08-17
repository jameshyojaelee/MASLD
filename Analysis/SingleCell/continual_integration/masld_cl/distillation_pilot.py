"""Control-only evaluation of the prospectively locked latent-distillation grid."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _bundle_centroids, _paired_improvement, _validate_harmony_info
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .control_evaluation import (
    _neighbor_jaccard_loss, _reference_f1_change_ci, _shift_and_standard_error,
)
from .embedding import load_embedding, matched_rows
from .firewall import assert_control_only_metric_names, validate_program_firewall
from .pilot import _run_identity


def load_distillation_policy(config: dict[str, Any], path: str | Path) -> dict[str, Any]:
    policy_path = Path(path).resolve()
    expected_path = Path(config["_config_path"]).resolve().parent / "reference" / policy_path.name
    if policy_path != expected_path or policy_path.name != "latent_distillation_policy_v3.json":
        raise ContractError("latent-distillation policy must be the source-controlled v3 policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    alphas = list(map(float, policy.get("alpha_grid", [])))
    if (
        policy.get("schema_version") != "masld-cl-latent-distillation-policy-v3"
        or policy.get("config_sha256") != config["_config_sha256"]
        or alphas != [100.0, 300.0, 1000.0, 3000.0, 10000.0]
        or float(policy.get("replay_fraction", -1)) != 0.20
        or int(policy.get("seed", -1)) != int(config["screen"]["seed"])
        or policy.get("case_stage_program_hero_gene_and_cas13_outcomes_locked") is not True
        or policy.get("zero_query_penalty") is not True
    ):
        raise ContractError("latent-distillation policy violates its prospective contract")
    parent_path = Path(config["_config_path"]).resolve().parent / policy[
        "parent_v2_decision_relative_path"
    ]
    if sha256_path(parent_path) != policy["parent_v2_decision_sha256"]:
        raise ContractError("v3 policy is not anchored to the immutable failed v2 result")
    with parent_path.open() as handle:
        parent = json.load(handle)
    if parent.get("passed") is not False or parent.get("outcomes_unlocked") is not False:
        raise ContractError("v3 parent was not a locked control-only failure")
    return policy


def _assert_policy_in_execution_source(record: dict[str, Any], policy_path: Path) -> None:
    with Path(record["execution_lock_realpath"]).open() as handle:
        lock = json.load(handle)
    relative = "reference/latent_distillation_policy_v3.json"
    matches = [
        entry for entry in lock["source_identity"]["entries"]
        if entry.get("path") == relative
    ]
    if len(matches) != 1 or matches[0].get("sha256") != sha256_path(policy_path):
        raise ContractError("execution source identity does not own the v3 policy")


def _candidate_metrics(
    config: dict[str, Any], reference, candidate, harmony_centroids,
    architecture_centroids, *, seed_offset: int,
) -> tuple[dict[str, float], dict[str, Any]]:
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
    centroids = _bundle_centroids(candidate)
    harmony = _paired_improvement(
        centroids, harmony_centroids, config["bootstrap"]["replicates"],
        config["bootstrap"]["seed"] + 2701,
    )
    architecture = _paired_improvement(
        centroids, architecture_centroids, config["bootstrap"]["replicates"],
        config["bootstrap"]["seed"] + 2801,
    )
    shift, shift_se = _shift_and_standard_error(
        candidate[2], candidate[1], config["bootstrap"]["replicates"],
        config["bootstrap"]["seed"] + 2901 + seed_offset,
    )
    metrics = {
        "reference_macro_f1_change_ci_low": f1["ci_low"],
        "reference_neighborhood_jaccard_loss": jaccard,
        "alignment_improvement_vs_harmony": harmony["improvement"],
        "alignment_improvement_vs_architecture_surgery": architecture["improvement"],
        "shift_control": shift,
        "shift_control_standard_error": shift_se,
    }
    assert_control_only_metric_names(list(metrics))
    return metrics, {"harmony": harmony, "architecture_surgery": architecture}


def evaluate_distillation_grid(
    config: dict[str, Any], policy_value: str | Path,
    reference_embedding: str | Path, reference_record: str | Path,
    harmony_embedding: str | Path,
    architecture_embedding: str | Path, architecture_record: str | Path,
    candidates: list[tuple[float, str, str]], output: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path = Path(policy_value).resolve()
    policy = load_distillation_policy(config, policy_path)
    expected_alphas = list(map(float, policy["alpha_grid"]))
    supplied = [float(value[0]) for value in candidates]
    if len(supplied) != len(set(supplied)) or set(supplied) != set(expected_alphas):
        raise ContractError("evaluation requires the exact prospectively locked alpha grid")

    reference, _, reference_execution, reference_run_path = _run_identity(
        config, reference_embedding, reference_record, "masld-cl-reference-v1"
    )
    architecture, architecture_run, architecture_execution, architecture_run_path = _run_identity(
        config, architecture_embedding, architecture_record, "masld-cl-update-v1"
    )
    expected_query = set(config["evaluation"]["powered_query_studies"])
    if (
        architecture_run.get("method") != "architecture_surgery"
        or architecture_run.get("production") is not False
        or architecture_run.get("seed") != policy["seed"]
        or float(architecture_run.get("ewc_lambda")) != 0
        or float(architecture_run.get("replay_fraction")) != 0
        or set(architecture_run.get("query_datasets", [])) != expected_query
    ):
        raise ContractError("architecture comparator has the wrong v3 context")
    harmony = load_embedding(harmony_embedding)
    _validate_harmony_info(harmony[0], config)
    source_ids = {
        reference_execution["source_identity_sha256"],
        architecture_execution["source_identity_sha256"],
    }
    _assert_policy_in_execution_source(reference_execution, policy_path)
    _assert_policy_in_execution_source(architecture_execution, policy_path)

    rows = []
    candidate_sources = []
    for index, (alpha, embedding_value, record_value) in enumerate(
        sorted(candidates, key=lambda value: value[0])
    ):
        candidate, run, execution, run_path = _run_identity(
            config, embedding_value, record_value, "masld-cl-update-v1"
        )
        distillation = run.get("distillation") or {}
        if (
            run.get("method") != "latent_distillation"
            or run.get("production") is not False
            or run.get("seed") != policy["seed"]
            or float(run.get("ewc_lambda")) != 0
            or float(run.get("replay_fraction")) != float(policy["replay_fraction"])
            or float(run.get("distillation_weight")) != float(alpha)
            or run.get("replay_mode") != "random"
            or set(run.get("query_datasets", [])) != expected_query
            or run.get("held_out_datasets")
            or run.get("fisher") is not None
            or distillation.get("target_statistic") != "qz.loc"
            or distillation.get("encoder_mode") != "evaluation"
            or distillation.get("anchor_after_reference_load_and_batch_expansion") is not True
            or int(distillation.get("n_targets", -1)) != int(run.get("n_replay_cells", -2))
            or int(distillation.get("n_query_targets", -1)) != 0
            or float(distillation.get("zero_at_anchor_maximum", 1)) > 1e-7
        ):
            raise ContractError(f"candidate violates v3 contract: alpha={alpha:g}")
        _assert_policy_in_execution_source(execution, policy_path)
        source_ids.add(execution["source_identity_sha256"])
        metrics, bootstrap = _candidate_metrics(
            config, reference, candidate, _bundle_centroids(harmony),
            _bundle_centroids(architecture), seed_offset=index,
        )
        threshold = config["gates"]["minimum_control_alignment_improvement"]
        gates = {
            "reference_macro_f1": metrics["reference_macro_f1_change_ci_low"]
            > -config["gates"]["reference_macro_f1_margin"],
            "reference_neighborhood": metrics["reference_neighborhood_jaccard_loss"]
            <= config["gates"]["reference_neighborhood_jaccard_loss"],
            "harmony_improvement": metrics["alignment_improvement_vs_harmony"] >= threshold
            and bootstrap["harmony"]["ci_low"] > 0,
            "architecture_improvement": metrics[
                "alignment_improvement_vs_architecture_surgery"
            ] >= threshold and bootstrap["architecture_surgery"]["ci_low"] > 0,
        }
        rows.append({
            "distillation_weight": float(alpha),
            "replay_fraction": float(policy["replay_fraction"]),
            "metrics": metrics, "bootstrap": bootstrap, "gates": gates,
            "eligible": all(gates.values()),
        })
        candidate_sources.append({
            "distillation_weight": float(alpha),
            "run_manifest": str(run_path), "run_manifest_sha256": sha256_path(run_path),
            "execution_record": str(Path(record_value).resolve()),
            "execution_record_sha256": sha256_path(record_value),
        })
    if len(source_ids) != 1:
        raise ContractError("v3 reference, comparator, and candidates are not source matched")

    eligible = [row for row in rows if row["eligible"]]
    selected = None
    if eligible:
        best = min(eligible, key=lambda row: row["metrics"]["shift_control"])
        within_one_se = best["metrics"]["shift_control"] + best["metrics"][
            "shift_control_standard_error"
        ]
        selected = min(
            row["distillation_weight"] for row in eligible
            if row["metrics"]["shift_control"] <= within_one_se
        )
    result = {
        "schema_version": "masld-cl-latent-distillation-pilot-v3",
        "config_sha256": config["_config_sha256"],
        "policy_realpath": str(policy_path), "policy_sha256": sha256_path(policy_path),
        "control_only": True, "complete_grid": True, "settings": rows,
        "eligible_weights": [row["distillation_weight"] for row in eligible],
        "selected_pilot_weight": selected,
        "passed": bool(eligible), "stop_required": not bool(eligible),
        "outcomes_unlocked": False,
        "source_identity_sha256": next(iter(source_ids)),
        "sources": {
            "reference_run": str(reference_run_path),
            "architecture_run": str(architecture_run_path),
            "harmony_embedding": str(Path(harmony_embedding).resolve()),
            "candidate_runs": candidate_sources,
        },
    }
    write_json_exclusive(output, result)
    return result
