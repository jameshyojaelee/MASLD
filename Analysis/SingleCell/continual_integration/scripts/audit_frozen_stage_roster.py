#!/usr/bin/env python3
"""Prove that correcting Liver Atlas donors cannot alter frozen stage scores."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from masld_cl.config import load_config, repo_path, write_json_exclusive
from masld_cl.contracts import ContractError
from masld_cl.firewall import sha256_file, validate_program_firewall


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    hashes = validate_program_firewall(config)
    relative = next(
        path for path in config["program_firewall"]
        if path.endswith("donor_program_scores_primary.tsv")
    )
    path = repo_path(config, relative)
    donors = defaultdict(set)
    rows = 0
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            donors[row["dataset"]].add(row["donor"])
    observed = {dataset: len(values) for dataset, values in sorted(donors.items())}
    expected = {"GSE174748": 4, "GSE185477": 3, "GSE189600": 2, "GSE202379": 37, "GSE244832": 18}
    if observed != expected or sum(observed.values()) != 64 or "Liver_Atlas" in observed:
        raise ContractError(f"frozen stage roster changed: {observed}")
    result = {
        "schema_version": "masld-cl-frozen-stage-roster-v1",
        "score_file": str(path),
        "score_file_sha256": sha256_file(path),
        "program_firewall_sha256": hashes,
        "rows": rows,
        "biological_donors": 64,
        "donors_by_dataset": observed,
        "liver_atlas_rows": 0,
        "stage_analysis_unchanged_by_liver_donor_correction": True,
    }
    write_json_exclusive(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
