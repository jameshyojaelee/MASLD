#!/usr/bin/env python3
"""Evaluate V17 after validating its reference-only selection lock."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from masld_cl.benchmark import _label_scores
from masld_cl.config import load_config, write_json_exclusive
from masld_cl.contracts import ContractError, sha256_path
from masld_cl.embedding import load_embedding


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--v14-outcomes", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    candidate_path = Path(args.candidate).resolve()
    candidate = load_embedding(candidate_path)
    manifest_path = candidate_path.parent / "reference_selected_dual_head_manifest.json"
    lock_path = candidate_path.parent / "reference_selection_lock.json"
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    with lock_path.open() as handle:
        lock = json.load(handle)
    if (
        manifest.get("schema_version") != "masld-cl-reference-selected-dual-head-v17"
        or lock.get("schema_version") != "masld-cl-reference-only-head-selection-lock-v17"
        or manifest.get("query_labels_used_for_selection") is not False
        or lock.get("query_audit_labels_read") is not False
        or manifest.get("selection_lock", {}).get("sha256") != sha256_path(lock_path)
        or manifest.get("weights", {}).get("adapted") != lock.get("selected_adapted_weight")
        or manifest.get("latent_bitwise_identical_to_v14") is not True
        or manifest.get("nonclassifier_digest_before") != manifest.get("nonclassifier_digest_after")
        or manifest.get("embedding") != candidate[0]
    ):
        raise ContractError("V17 manifest or reference-only selection lock differs")
    baseline_path = Path(args.baseline).resolve()
    with baseline_path.open() as handle:
        authority = json.load(handle)
    best = authority["best_non_cl_method"]
    baseline = authority["method_scores"][best]
    overall, lineages = _label_scores(candidate, config["lineages"])
    overall_change = overall - baseline["mean_donor_macro_f1"]
    lineage_changes = {
        lineage: lineages[lineage] - baseline["mean_donor_lineage_f1"][lineage]
        for lineage in config["lineages"]
    }
    label_pass = bool(
        overall_change >= -config["gates"]["query_macro_f1_margin"]
        and min(lineage_changes.values()) >= -config["gates"]["major_lineage_f1_margin"]
    )
    v14_path = Path(args.v14_outcomes).resolve()
    with v14_path.open() as handle:
        v14 = json.load(handle)
    if v14.get("disease_preservation", {}).get("pass") is not True:
        raise ContractError("V17 cannot inherit a failed disease gate")
    result = {
        "schema_version": "masld-cl-reference-selected-dual-head-label-diagnostic-v17",
        "config_sha256": config["_config_sha256"],
        "candidate": {"path": str(candidate_path), "sha256": sha256_path(candidate_path)},
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256_path(manifest_path)},
        "selection_lock": {"path": str(lock_path.resolve()), "sha256": sha256_path(lock_path)},
        "baseline": {"path": str(baseline_path), "sha256": sha256_path(baseline_path)},
        "best_non_cl_method": best,
        "selected_adapted_weight": lock["selected_adapted_weight"],
        "candidate_mean_donor_macro_f1": overall,
        "baseline_mean_donor_macro_f1": baseline["mean_donor_macro_f1"],
        "query_macro_f1_change": overall_change,
        "candidate_mean_donor_lineage_f1": lineages,
        "major_lineage_f1_changes": lineage_changes,
        "query_label_concordance_pass": label_pass,
        "disease_preservation_inherited": {
            "path": str(v14_path), "sha256": sha256_path(v14_path), "pass": True
        },
        "pilot_pass": label_pass,
        "development_only_until_independent_held_study_confirmation": True,
        "labels_are_audit_targets_and_were_not_overwritten": True,
    }
    write_json_exclusive(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
