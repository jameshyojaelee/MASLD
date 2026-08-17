#!/usr/bin/env python3
"""Re-derive the frozen 117-program integrity contract."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from masld_cl.config import load_config, repo_path, write_json_exclusive
from masld_cl.contracts import ContractError, sha256_path
from masld_cl.firewall import validate_program_firewall


def _rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def audit(config):
    observed = validate_program_firewall(config)
    files = {Path(path).name: repo_path(config, path) for path in observed}
    registry = _rows(files["program_registry_v2.tsv"])
    membership = _rows(files["program_membership_v2.tsv"])
    primary = _rows(files["donor_program_scores_primary.tsv"])
    equal_run = _rows(files["donor_program_scores_equal_run.tsv"])
    specification = _rows(files["analysis_specification.tsv"])
    if len(registry) != 117 or len({row["program_uid"] for row in registry}) != 117:
        raise ContractError("frozen program registry does not contain 117 unique programs")
    if len({(row["cell_type"], row["module"]) for row in registry}) != 117:
        raise ContractError("frozen program cell-type/module keys are not unique")
    registry_by_uid = {row["program_uid"]: row for row in registry}
    members = defaultdict(list)
    for row in membership:
        if row["program_uid"] not in registry_by_uid:
            raise ContractError("membership contains a program absent from the registry")
        members[row["program_uid"]].append(row)
    if set(members) != set(registry_by_uid):
        raise ContractError("membership lacks a frozen registry program")
    for uid, rows in members.items():
        registry_row = registry_by_uid[uid]
        if (
            len(rows) != int(registry_row["n_source_genes"])
            or {row["membership_sha256"] for row in rows}
            != {registry_row["membership_sha256"]}
            or any(float(row["source_weight"]) <= 0 for row in rows)
            or abs(sum(float(row["original_l1_weight"]) for row in rows) - 1.0) > 1e-9
        ):
            raise ContractError(f"frozen membership or weights differ: {uid}")
    score_summaries = {}
    for name, rows in (("primary", primary), ("equal_run", equal_run)):
        if {row["program_uid"] for row in rows} != set(registry_by_uid):
            raise ContractError(f"{name} donor scores lack a frozen program")
        per_program = Counter(row["program_uid"] for row in rows)
        score_summaries[name] = {
            "rows": len(rows),
            "programs": len(per_program),
            "minimum_donor_scores_per_program": min(per_program.values()),
            "maximum_donor_scores_per_program": max(per_program.values()),
            "score_file_sha256": sha256_path(files[
                "donor_program_scores_primary.tsv"
                if name == "primary" else "donor_program_scores_equal_run.tsv"
            ]),
        }
    contract = config["program_inference_contract"]
    if (
        len(specification) != 1
        or specification[0]["primary_model"] != contract["primary_model"]
        or specification[0]["external_outcomes_read"] != "FALSE"
        or contract["robust_variance"] != "HC3"
        or contract["complete_program_registry_size"] != 117
        or contract["weighted_bh_family_size"] != 113
        or contract["unweighted_bh_family_size"] != 114
        or contract["new_discovery_allowed"] is not False
    ):
        raise ContractError("frozen program inference specification differs")
    return {
        "schema_version": "masld-cl-program-firewall-audit-v39",
        "config_sha256": config["_config_sha256"],
        "passed": True,
        "observed_file_sha256": observed,
        "registry": {
            "programs": len(registry),
            "unique_program_uids": len(registry_by_uid),
            "membership_rows": len(membership),
            "membership_file_sha256": sha256_path(files["program_membership_v2.tsv"]),
            "names_weights_and_memberships_byte_frozen": True,
        },
        "scores": score_summaries,
        "inference": {
            "primary_model": contract["primary_model"],
            "robust_variance": contract["robust_variance"],
            "weighted_bh_family_size": contract["weighted_bh_family_size"],
            "unweighted_bh_family_size": contract["unweighted_bh_family_size"],
            "new_discovery_allowed": contract["new_discovery_allowed"],
            "external_outcomes_read": False,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = audit(load_config(args.config))
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "program_firewall_audit.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
