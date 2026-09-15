#!/usr/bin/env python3
"""Bind one job-local gzip preflight to its completed GSE105127 v2 stage."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


STAGE_STATUS = {
    "rsem_reference": "passed",
    "rrbs_collapse": "complete",
    "rna_quantification": "complete",
    "cpg_crosswalk": "passed_failure_aware_roundtrip",
}


class GSE105127StageLinkV2Error(ValueError):
    """Raised when the preflight and completed stage cannot be bound."""


def plan_rows(plan_root: Path, stage: str, bundle_id: int) -> list[dict[str, str]]:
    filename = "rrbs_rows.tsv" if stage == "rrbs_collapse" else "rna_rows.tsv"
    with (plan_root / filename).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        forbidden = {"phenotype", "label", "outcome", "disease", "fibrosis", "nas", "sex", "age", "bmi"}
        if forbidden & set(reader.fieldnames or ()):
            raise GSE105127StageLinkV2Error("source plan violates the outcome firewall")
        rows = [row for row in reader if int(row["bundle_id"]) == bundle_id]
    if not rows:
        raise GSE105127StageLinkV2Error("stage-link bundle is empty")
    return rows


def freeze_link(
    *,
    stage: str,
    plan_root: Path,
    preflight_root: Path,
    stage_root: Path,
    output: Path,
    bundle_id: int | None,
) -> dict[str, object]:
    if output.exists():
        raise GSE105127StageLinkV2Error("refusing to overwrite stage link")
    if stage not in STAGE_STATUS:
        raise GSE105127StageLinkV2Error("unsupported stage link")
    verify_frozen_tree(plan_root)
    preflight_manifest = verify_frozen_tree(preflight_root)
    preflight = json.loads((preflight_root / "receipt.json").read_text(encoding="utf-8"))
    slurm_job_id = os.environ.get("SLURM_JOB_ID")
    slurm_array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
    if (
        slurm_job_id is None
        or preflight.get("status") != "passed_exactly_one_gzip_member_per_input"
        or preflight.get("stage") != stage
        or preflight.get("bundle_id") != bundle_id
        or preflight.get("slurm_job_id") != slurm_job_id
        or preflight.get("slurm_array_task_id") != slurm_array_task_id
        or preflight.get("all_gzip_test_passed") is not True
        or preflight.get("all_exactly_one_member") is not True
        or preflight.get("molecular_values_parsed") is not False
        or preflight.get("labels_accessed") is not False
        or preflight.get("fit_or_score_performed") is not False
    ):
        raise GSE105127StageLinkV2Error("job-local preflight receipt differs")
    logical_stage = stage_root
    if stage in {"rrbs_collapse", "rna_quantification"}:
        if bundle_id is None:
            raise GSE105127StageLinkV2Error("bundle stage requires a bundle ID")
        logical_stage = stage_root / f"bundle_{bundle_id:02d}"
    stage_manifest = verify_frozen_tree(logical_stage)
    stage_receipt = json.loads((logical_stage / "receipt.json").read_text(encoding="utf-8"))
    if (
        stage_receipt.get("status") != STAGE_STATUS[stage]
        or stage_receipt.get("labels_accessed") is not False
        or stage_receipt.get("fit_or_score_performed") is not False
        or (bundle_id is not None and stage_receipt.get("bundle_id") != bundle_id)
    ):
        raise GSE105127StageLinkV2Error("completed stage receipt differs")
    member_artifacts: dict[str, str] = {}
    if stage in {"rrbs_collapse", "rna_quantification"}:
        assert bundle_id is not None
        for row in plan_rows(plan_root, stage, bundle_id):
            participant = stage_root / "participants" / row["row_id"]
            verify_frozen_tree(participant)
            receipt = json.loads((participant / "receipt.json").read_text(encoding="utf-8"))
            expected_status = "passed" if stage == "rrbs_collapse" else "passed_raw_scale"
            if (
                receipt.get("status") != expected_status
                or receipt.get("row_id") != row["row_id"]
                or receipt.get("participant_group_id") != row["participant_group_id"]
                or receipt.get("zone") != row["zone"]
                or receipt.get("labels_accessed") is not False
                or receipt.get("fit_or_score_performed") is not False
            ):
                raise GSE105127StageLinkV2Error("completed participant receipt differs")
            member_artifacts[row["row_id"]] = sha256_file(participant / "ARTIFACTS.json")
    receipt = {
        "schema_version": "masld-bench-gse105127-gzip-stage-link-v2",
        "status": "job_local_preflight_bound_to_completed_stage",
        "stage": stage,
        "bundle_id": bundle_id,
        "slurm_job_id": slurm_job_id,
        "slurm_array_task_id": slurm_array_task_id,
        "preflight_artifacts_sha256": sha256_file(preflight_root / "ARTIFACTS.json"),
        "preflight_file_count": preflight["file_count"],
        "stage_artifacts_sha256": sha256_file(logical_stage / "ARTIFACTS.json"),
        "member_artifacts_sha256": member_artifacts,
        "molecular_values_opened_by_linker": False,
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
    }
    output.mkdir(mode=0o750)
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse105127_gzip_stage_link_v2",
            "stage": stage,
            "bundle_id": bundle_id,
            "slurm_job_id": slurm_job_id,
            "preflight_manifest_sha256": sha256_file(preflight_root / "ARTIFACTS.json"),
            "stage_manifest_sha256": sha256_file(logical_stage / "ARTIFACTS.json"),
            "status": "passed",
        },
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=tuple(STAGE_STATUS))
    parser.add_argument("--plan-root", required=True, type=Path)
    parser.add_argument("--preflight-root", required=True, type=Path)
    parser.add_argument("--stage-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bundle-id", type=int)
    args = parser.parse_args()
    print(json.dumps(freeze_link(**vars(args)), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
