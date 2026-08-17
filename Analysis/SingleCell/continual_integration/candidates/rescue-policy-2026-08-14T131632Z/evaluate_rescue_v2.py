#!/usr/bin/env python3
"""Evaluate the prospectively locked, control-only stronger-lambda rescue."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from masld_cl.benchmark import _bundle_centroids, _paired_improvement, _validate_harmony_info
from masld_cl.config import load_config, write_json_exclusive
from masld_cl.contracts import ContractError, sha256_path
from masld_cl.control_evaluation import (
    _neighbor_jaccard_loss,
    _reference_f1_change_ci,
    _shift_and_standard_error,
)
from masld_cl.embedding import load_embedding, matched_rows
from masld_cl.execution import require_execution_ownership, verify_execution_record
from masld_cl.firewall import assert_control_only_metric_names, validate_program_firewall


def _load_policy(path: Path, config: dict) -> dict:
    with path.open() as handle:
        policy = json.load(handle)
    parent = policy.get("parent", {})
    failed = Path(parent.get("failed_pilot_realpath", "")).resolve()
    if (
        policy.get("schema_version") != "masld-cl-control-only-rescue-policy-v2"
        or policy.get("status") != "locked_before_v2_training"
        or policy.get("decision", {}).get("outcomes_unlocked") is not False
        or Path(parent.get("config_realpath", "")).resolve()
        != Path(config["_config_path"]).resolve()
        or sha256_path(config["_config_path"]) != parent.get("config_sha256")
        or sha256_path(failed) != parent.get("failed_pilot_sha256")
    ):
        raise ContractError("rescue policy or its failed-v1 parent is invalid")
    with failed.open() as handle:
        parent_result = json.load(handle)
    if parent_result.get("passed") is not False or parent_result.get("control_only") is not True:
        raise ContractError("rescue policy is not anchored to a failed control-only v1 pilot")
    if policy["rationale"].get("objective_change") != "none":
        raise ContractError("rescue may change only lambda, not the EWC objective")
    return policy


def _load_owned_run(config: dict, embedding_value: str, record_value: str, schema: str):
    embedding_path = Path(embedding_value).resolve()
    bundle = load_embedding(embedding_path)
    manifest_name = "reference_manifest.json" if schema == "masld-cl-reference-v1" else "update_manifest.json"
    manifest_path = embedding_path.parent / manifest_name
    with manifest_path.open() as handle:
        run = json.load(handle)
    if (
        run.get("schema_version") != schema
        or run.get("config_sha256") != config["_config_sha256"]
        or run.get("model_kind") != "all_lineage"
        or run.get("embedding") != bundle[0]
    ):
        raise ContractError(f"training manifest does not own embedding: {embedding_path}")
    record_path = Path(record_value).resolve()
    record = verify_execution_record(
        record_path, Path(config["_config_path"]).resolve().parent,
        config["_config_sha256"],
    )
    require_execution_ownership(record, [embedding_path, manifest_path], role="rescue pilot")
    return bundle, run, record, manifest_path


def _parse_candidate(value: str) -> tuple[float, str, str]:
    fields = value.split("|", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("candidate must be LAMBDA|EMBEDDING|EXECUTION_RECORD")
    return float(fields[0]), fields[1], fields[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--reference-execution-record", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--architecture-embedding", required=True)
    parser.add_argument("--architecture-execution-record", required=True)
    parser.add_argument("--candidate", action="append", type=_parse_candidate, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    validate_program_firewall(config)
    policy_path = Path(args.policy).resolve()
    policy = _load_policy(policy_path, config)
    expected_lambdas = list(map(float, policy["pilot"]["lambda_values"]))
    supplied_lambdas = [value[0] for value in args.candidate]
    if len(supplied_lambdas) != len(set(supplied_lambdas)) or set(supplied_lambdas) != set(expected_lambdas):
        raise ContractError("rescue evaluation requires the exact prospectively locked lambda grid")

    reference, _, reference_record, reference_run_path = _load_owned_run(
        config, args.reference_embedding, args.reference_execution_record,
        "masld-cl-reference-v1",
    )
    architecture, architecture_run, architecture_record, architecture_run_path = _load_owned_run(
        config, args.architecture_embedding, args.architecture_execution_record,
        "masld-cl-update-v1",
    )
    if (
        architecture_run.get("method") != "architecture_surgery"
        or architecture_run.get("production") is not False
        or architecture_run.get("seed") != policy["pilot"]["seed"]
        or float(architecture_run.get("ewc_lambda")) != 0
        or float(architecture_run.get("replay_fraction")) != 0
    ):
        raise ContractError("architecture comparator has the wrong context")
    harmony = load_embedding(args.harmony_embedding)
    _validate_harmony_info(harmony[0], config)
    reference_cells = reference[2].reset_index(drop=True)
    harmony_centroids = _bundle_centroids(harmony)
    architecture_centroids = _bundle_centroids(architecture)
    policy_sha = sha256_path(policy_path)
    gates = policy["gates"]
    rows = []
    sources = []

    for index, (ewc_lambda, embedding_value, record_value) in enumerate(
        sorted(args.candidate, key=lambda value: value[0])
    ):
        candidate, run, record, run_path = _load_owned_run(
            config, embedding_value, record_value, "masld-cl-update-v1"
        )
        fisher = run.get("fisher") or {}
        spec_path = Path(record["spec_realpath"])
        with spec_path.open() as handle:
            spec = json.load(handle)
        if (
            run.get("method") != "continual_learning"
            or run.get("replay_mode") != "random"
            or run.get("sensitivity_only") is not False
            or run.get("production") is not False
            or run.get("seed") != policy["pilot"]["seed"]
            or float(run.get("ewc_lambda")) != ewc_lambda
            or float(run.get("replay_fraction")) != float(policy["pilot"]["replay_fraction"])
            or set(run.get("query_datasets", [])) != set(policy["pilot"]["query_datasets"])
            or run.get("held_out_datasets")
            or set(fisher.get("control_fisher_datasets", []))
            != set(policy["pilot"]["control_fisher_datasets"])
            or spec.get("rescue_policy_realpath") != str(policy_path)
            or spec.get("rescue_policy_sha256") != policy_sha
        ):
            raise ContractError(f"candidate does not match rescue policy: lambda={ewc_lambda:g}")

        candidate_reference = candidate[2].loc[candidate[2]["strict_reference"]].reset_index(drop=True)
        left, right = matched_rows(reference_cells, candidate_reference)
        f1 = _reference_f1_change_ci(
            reference_cells.iloc[left].reset_index(drop=True),
            candidate_reference.iloc[right].reset_index(drop=True),
            policy["pilot"]["bootstrap_replicates"], config["bootstrap"]["seed"],
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
        harmony_result = _paired_improvement(
            centroids, harmony_centroids, policy["pilot"]["bootstrap_replicates"],
            config["bootstrap"]["seed"] + 1701,
        )
        architecture_result = _paired_improvement(
            centroids, architecture_centroids, policy["pilot"]["bootstrap_replicates"],
            config["bootstrap"]["seed"] + 1801,
        )
        shift, shift_se = _shift_and_standard_error(
            candidate[2], candidate[1], policy["pilot"]["bootstrap_replicates"],
            config["bootstrap"]["seed"] + 1901 + index,
        )
        metrics = {
            "reference_macro_f1_change_ci_low": f1["ci_low"],
            "reference_neighborhood_jaccard_loss": jaccard,
            "alignment_improvement_vs_harmony": harmony_result["improvement"],
            "alignment_improvement_vs_architecture_surgery": architecture_result["improvement"],
            "shift_control": shift,
            "shift_control_standard_error": shift_se,
        }
        assert_control_only_metric_names(list(metrics))
        setting_gates = {
            "reference_macro_f1": metrics["reference_macro_f1_change_ci_low"]
            > gates["reference_macro_f1_change_ci_low_greater_than"],
            "reference_neighborhood": metrics["reference_neighborhood_jaccard_loss"]
            <= gates["reference_neighborhood_jaccard_loss_at_most"],
            "harmony_noncatastrophic": metrics["alignment_improvement_vs_harmony"]
            >= gates["alignment_improvement_vs_harmony_at_least"],
            "architecture_noncatastrophic": metrics["alignment_improvement_vs_architecture_surgery"]
            >= gates["alignment_improvement_vs_architecture_surgery_at_least"],
        }
        rows.append({
            "ewc_lambda": ewc_lambda,
            "replay_fraction": policy["pilot"]["replay_fraction"],
            "metrics": metrics,
            "bootstrap": {"harmony": harmony_result, "architecture_surgery": architecture_result},
            "gates": setting_gates,
            "eligible": all(setting_gates.values()),
        })
        sources.append({
            "ewc_lambda": ewc_lambda,
            "run_manifest": str(run_path),
            "run_manifest_sha256": sha256_path(run_path),
            "execution_record": str(Path(record_value).resolve()),
            "execution_record_sha256": sha256_path(record_value),
        })

    eligible = [row for row in rows if row["eligible"]]
    selected = None
    if eligible:
        best = min(eligible, key=lambda row: row["metrics"]["shift_control"])
        threshold = best["metrics"]["shift_control"] + best["metrics"]["shift_control_standard_error"]
        selected = min(
            (row for row in eligible if row["metrics"]["shift_control"] <= threshold),
            key=lambda row: row["ewc_lambda"],
        )["ewc_lambda"]
    result = {
        "schema_version": "masld-cl-control-only-rescue-result-v2",
        "config_sha256": config["_config_sha256"],
        "policy_realpath": str(policy_path),
        "policy_sha256": policy_sha,
        "control_only": True,
        "complete_grid": True,
        "settings": rows,
        "eligible_lambdas": [row["ewc_lambda"] for row in eligible],
        "selected_pilot_scale": selected,
        "passed": bool(eligible),
        "stop_required": not bool(eligible),
        "outcomes_unlocked": False,
        "sources": {
            "reference_run": {
                "path": str(reference_run_path), "sha256": sha256_path(reference_run_path),
                "execution_record": str(Path(args.reference_execution_record).resolve()),
                "execution_record_sha256": sha256_path(args.reference_execution_record),
                "output_tree_sha256": reference_record["result_identity"]["output_tree_sha256"],
            },
            "architecture_run": {
                "path": str(architecture_run_path), "sha256": sha256_path(architecture_run_path),
                "execution_record": str(Path(args.architecture_execution_record).resolve()),
                "execution_record_sha256": sha256_path(args.architecture_execution_record),
                "output_tree_sha256": architecture_record["result_identity"]["output_tree_sha256"],
            },
            "harmony_embedding": {"path": str(Path(args.harmony_embedding).resolve()), "sha256": sha256_path(args.harmony_embedding)},
            "candidate_runs": sources,
            "evaluator": {"path": str(Path(__file__).resolve()), "sha256": sha256_path(__file__)},
        },
    }
    write_json_exclusive(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
