"""Locked best-non-CL label-concordance benchmark for the V9 candidate."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .benchmark import _label_scores
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding, matched_rows
from .execution import require_execution_ownership, verify_execution_record
from .firewall import load_retargeted_bridge_selection_lock, validate_program_firewall


METHODS = {"architecture_surgery", "de_novo", "fine_tune", "replay_only", "ewc_only"}


def _mapping(values: list[str], expected: set[str], name: str) -> dict[str, Path]:
    result = {}
    for value in values:
        key, separator, path = value.partition("=")
        if separator != "=" or key in result:
            raise ContractError(f"invalid or duplicate {name} mapping: {value}")
        result[key] = Path(path).resolve()
    if set(result) != expected:
        raise ContractError(f"{name} mappings require exactly {sorted(expected)}")
    return result


def evaluate_best_non_cl_labels(
    config: dict[str, Any], policy_value: str | Path, selection_lock: str | Path,
    candidate_embedding: str | Path, embedding_values: list[str],
    execution_record_values: list[str], output: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    policy_path = Path(policy_value).resolve()
    expected_policy = pipeline_root / "reference" / "non_cl_baseline_policy_v10.json"
    if policy_path != expected_policy:
        raise ContractError("non-CL benchmark requires the source-controlled V10 policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    selection_path = Path(selection_lock).resolve()
    selection = load_retargeted_bridge_selection_lock(selection_path, config)
    if (
        policy.get("schema_version") != "masld-cl-non-cl-baseline-policy-v10"
        or policy.get("config_sha256") != config["_config_sha256"]
        or sha256_path(selection_path) != policy.get("selection_lock_file_sha256")
        or selection.get("lock_sha256") != policy.get("selection_lock_sha256")
    ):
        raise ContractError("non-CL baseline policy or selection lock differs")
    embeddings = _mapping(embedding_values, METHODS, "embedding")
    records = _mapping(execution_record_values, METHODS, "execution-record")
    policy_settings = {
        row["method"]: (float(row["ewc_lambda"]), float(row["replay_fraction"]))
        for row in policy.get("baselines", [])
    }
    if set(policy_settings) != METHODS - {"architecture_surgery"}:
        raise ContractError("non-CL policy method roster is incomplete")
    candidate_path = Path(candidate_embedding).resolve()
    selected_all = None
    selected_architecture = None
    for source in selection["control_results"]:
        with Path(source["path"]).open() as handle:
            row = json.load(handle)
        if (
            row["model_kind"] == "all_lineage"
            and float(row["global_reference_weight"]) == float(selection["selected_global_reference_weight"])
        ):
            selected_all = Path(row["sources"]["candidate_embedding"]).resolve()
            selected_architecture = Path(row["sources"]["architecture_embedding"]).resolve()
    if selected_all != candidate_path:
        raise ContractError("label benchmark candidate is not the locked all-lineage V9 embedding")
    candidate = load_embedding(candidate_path)
    method_scores, sources = {}, []
    pipeline_root = Path(config["_config_path"]).resolve().parent
    for method in sorted(METHODS):
        path = embeddings[method]
        bundle = load_embedding(path)
        candidate_rows, baseline_rows = matched_rows(candidate[2], bundle[2])
        if len(candidate_rows) != len(candidate[2]) or len(baseline_rows) != len(bundle[2]):
            raise ContractError(f"non-CL comparator has a different cell roster: {method}")
        manifest_path = path.parent / ("denovo_manifest.json" if method == "de_novo" else "update_manifest.json")
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        expected_schema = "masld-cl-denovo-v1" if method == "de_novo" else "masld-cl-update-v1"
        if (
            manifest.get("schema_version") != expected_schema
            or manifest.get("config_sha256") != config["_config_sha256"]
            or manifest.get("model_kind") != "all_lineage"
            or manifest.get("seed") != int(policy["seed"])
            or set(manifest.get("query_datasets", [])) != set(policy["query_datasets"])
            or manifest.get("embedding") != bundle[0]
        ):
            raise ContractError(f"non-CL comparator provenance differs: {method}")
        if method == "architecture_surgery":
            if path != selected_architecture:
                raise ContractError("architecture comparator is not the exact locked V9 source")
        elif manifest.get("selection_lock_sha256") != selection["lock_sha256"]:
            raise ContractError(f"non-CL comparator predates the V9 outcome lock: {method}")
        if method != "de_novo" and manifest.get("method") != method:
            raise ContractError(f"non-CL update method differs: {method}")
        if method == "architecture_surgery":
            expected_setting = (0.0, 0.0)
        elif method != "de_novo":
            expected_setting = policy_settings[method]
        else:
            expected_setting = None
        if expected_setting is not None and (
            float(manifest.get("ewc_lambda", -1)) != expected_setting[0]
            or float(manifest.get("replay_fraction", -1)) != expected_setting[1]
        ):
            raise ContractError(f"non-CL comparator setting differs from policy: {method}")
        record = verify_execution_record(
            records[method], pipeline_root, config["_config_sha256"], allow_historical_source=True
        )
        require_execution_ownership(record, [path, manifest_path], role=f"non-CL {method}")
        overall, lineages = _label_scores(bundle, config["lineages"])
        method_scores[method] = {"mean_donor_macro_f1": overall, "mean_donor_lineage_f1": lineages}
        sources.append({
            "method": method, "embedding": str(path), "embedding_sha256": sha256_path(path),
            "training_manifest": str(manifest_path.resolve()),
            "training_manifest_sha256": sha256_path(manifest_path),
            "execution_record": str(records[method]), "execution_record_sha256": sha256_path(records[method]),
        })
    best = sorted(method_scores, key=lambda key: (-method_scores[key]["mean_donor_macro_f1"], key))[0]
    candidate_overall, candidate_lineages = _label_scores(candidate, config["lineages"])
    overall_change = candidate_overall - method_scores[best]["mean_donor_macro_f1"]
    lineage_changes = {
        lineage: candidate_lineages[lineage] - method_scores[best]["mean_donor_lineage_f1"][lineage]
        for lineage in config["lineages"]
    }
    passed = bool(
        overall_change >= -config["gates"]["query_macro_f1_margin"]
        and min(lineage_changes.values()) >= -config["gates"]["major_lineage_f1_margin"]
    )
    result = {
        "schema_version": "masld-cl-best-non-cl-label-benchmark-v10",
        "config_sha256": config["_config_sha256"],
        "policy": str(policy_path), "policy_sha256": sha256_path(policy_path),
        "selection_lock": str(selection_path), "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"],
        "candidate_embedding": str(candidate_path), "candidate_embedding_sha256": sha256_path(candidate_path),
        "method_scores": method_scores, "best_non_cl_method": best,
        "candidate_mean_donor_macro_f1": candidate_overall,
        "query_macro_f1_change_vs_best_non_cl": overall_change,
        "major_lineage_f1_changes_vs_best_non_cl": lineage_changes,
        "thresholds": {
            "query_macro_f1_change": -config["gates"]["query_macro_f1_margin"],
            "major_lineage_f1_change": -config["gates"]["major_lineage_f1_margin"],
        },
        "query_label_concordance_pass": passed, "sources": sources,
        "labels_are_audit_targets_and_were_not_overwritten": True,
    }
    write_json_exclusive(output, result)
    return result
