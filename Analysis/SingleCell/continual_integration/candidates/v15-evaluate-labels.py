#!/usr/bin/env python3
"""Evaluate the reference-head diagnostic against the frozen non-CL authority."""

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
    head_path = candidate_path.parent / "reference_head_manifest.json"
    with head_path.open() as handle:
        head = json.load(handle)
    if (
        head.get("schema_version") != "masld-cl-reference-head-preservation-v15"
        or head.get("query_labels_used") is not False
        or head.get("latent_bitwise_identical_to_v14") is not True
        or head.get("nonclassifier_digest_before")
        != head.get("nonclassifier_digest_after")
        or head.get("classifier_source_digest")
        != head.get("classifier_digest_after")
        or head.get("embedding") != candidate[0]
    ):
        raise ContractError("V15 reference-head manifest differs")
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
        and min(lineage_changes.values())
        >= -config["gates"]["major_lineage_f1_margin"]
    )
    v14_path = Path(args.v14_outcomes).resolve()
    with v14_path.open() as handle:
        v14 = json.load(handle)
    if v14.get("disease_preservation", {}).get("pass") is not True:
        raise ContractError("V15 cannot inherit a failed V14 disease gate")
    result = {
        "schema_version": "masld-cl-reference-head-label-diagnostic-v15",
        "config_sha256": config["_config_sha256"],
        "candidate": {"path": str(candidate_path), "sha256": sha256_path(candidate_path)},
        "reference_head_manifest": {"path": str(head_path.resolve()), "sha256": sha256_path(head_path)},
        "baseline": {"path": str(baseline_path), "sha256": sha256_path(baseline_path)},
        "best_non_cl_method": best,
        "candidate_mean_donor_macro_f1": overall,
        "baseline_mean_donor_macro_f1": baseline["mean_donor_macro_f1"],
        "query_macro_f1_change": overall_change,
        "candidate_mean_donor_lineage_f1": lineages,
        "major_lineage_f1_changes": lineage_changes,
        "query_label_concordance_pass": label_pass,
        "disease_preservation_inherited": {
            "path": str(v14_path),
            "sha256": sha256_path(v14_path),
            "pass": True,
            "basis": "V15 latent file is bitwise identical to the selected V14 latent file"
        },
        "pilot_pass": label_pass,
        "labels_are_audit_targets_and_were_not_overwritten": True
    }
    write_json_exclusive(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
