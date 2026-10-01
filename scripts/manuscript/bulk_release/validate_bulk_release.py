#!/usr/bin/env python3
"""Independent structural and statistical checks for a bulk-only release."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

from build_bulk_release import (RELEASE_ID, load_expected, require, sha256, validate_f0, validate_pooled,
                                validate_stage)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=Path, required=True)
    parser.add_argument("--release-id", default=RELEASE_ID)
    args = parser.parse_args()
    release_id = args.release_id
    root = args.release.resolve()
    require(root.name == release_id, f"Unexpected bulk-release directory: {root}")
    expected = load_expected(root / "artifacts/pooled/expected_counts.json")
    manifest = read_rows(root / "release_manifest.tsv")
    require(len(manifest) == 1, "Release manifest must have exactly one row")
    row = manifest[0]
    require(row == {
        "release_id": release_id,
        "release_state": "adopted_bulk_only",
        "scope": "five_cohort_bulk_dependents",
        "pooled_samples": str(expected["n_samples"]),
        "controls": str(expected["n_control"]),
        "disease": str(expected["n_disease"]),
        "tested_genes": str(expected["n_genes"]),
        "resource_synchronized_release": "pending",
    }, "Bulk-release manifest contents drift")
    release_json = json.loads((root / "BULK_RELEASE.json").read_text())
    require(release_json["release_id"] == release_id and release_json["resource_synchronized_release"] == "pending",
            "Bulk release JSON contradicts the adoption boundary")
    require("blocked" in release_json["figure3F"], "Figure 3F dependency state is not explicit")

    inventory = read_rows(root / "artifact_inventory.tsv")
    require(inventory, "Bulk release has no artifact inventory")
    release_paths = set()
    for item in inventory:
        relative = Path(item["release_path"])
        require(not relative.is_absolute() and ".." not in relative.parts, f"Invalid release-relative path: {relative}")
        artifact = root / relative
        require(artifact.is_file(), f"Inventoried artifact is missing: {artifact}")
        require(sha256(artifact) == item["sha256"], f"Inventoried hash mismatch: {artifact}")
        require(item["release_path"] not in release_paths, f"Duplicate inventory row: {relative}")
        release_paths.add(item["release_path"])

    pooled = validate_pooled(root / "artifacts/pooled/deg_results.csv", expected)
    stage = validate_stage(root / "artifacts/stage/stage_extension_all_gene_results.tsv", expected["n_genes"])
    f0 = validate_f0(root / "artifacts/f0_sensitivity/f0_arm_results.tsv.gz",
                     root / "artifacts/f0_sensitivity/f0_arm_summary.tsv", expected)
    family_rows = read_rows(root / "validation/statistical_family_checks.tsv")
    n_families = stage["stage_families"] + f0["f0_families"]
    require(len(family_rows) == n_families,
            f"Expected {n_families} recorded statistical families, found {len(family_rows)}")
    require((root / "figures/figure3/fig3e_stage_remodeling.pdf").is_file(), "Missing Figure 3E")
    require((root / "figures/figureS3/figs3_stage_remodeling_full.pdf").is_file(), "Missing Figure S3 stage panel")
    require(not (root / "figures/figure3/fig3f_stage_deg_genetics_matrix.pdf").exists(),
            "Figure 3F must remain outside the bulk-only release")
    scope = subprocess.run(["python3", "scripts/manuscript/validate_resource_scope.py"],
                           cwd=root.parents[3], text=True, capture_output=True)
    require(scope.returncode == 0, f"Resource-scope check failed: {scope.stdout}{scope.stderr}")
    report = {
        "release_id": release_id,
        "status": "PASS",
        "inventoried_artifacts": len(inventory),
        "pooled": pooled,
        "stage_families": stage["stage_families"],
        "f0_estimable_families": f0["f0_families"],
        "resource_scope": "PASS",
    }
    report_path = root / "validation/independent_validation.json"
    require(not report_path.exists(), f"Refusing to overwrite independent validation: {report_path}")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
