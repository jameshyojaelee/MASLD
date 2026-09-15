#!/usr/bin/env python3
"""Build and freeze the five-seed Cobolt v10 GPU dispatch bundle.

This is the recovery revision of ``build_cobolt_v9_gpu_bundle.py``.  GPU job
21100294 executed all 25 v9 runs and failed every one of them at run inclusion,
before any adapter action, because concurrent Figure 4 Resource-authority
adoption changed protected files while the job held the frozen v9 separation
baseline (config/artifacts/incidents/cobolt_v9_resource_authority_firewall_21100294.json).
The frozen ``continuation_policy`` forbids mid-run resume and requires a new
campaign revision, a distinct read-only bundle, and a distinct queue item, so
this builder refuses anything but the v10 campaign and the v10 bundle identity.

Unlike the v9 builder this one does not emit the GPU wrapper.  The wrapper is a
reviewed, checked-in file at slurm/run_cobolt_5seed_25run_v10_generic7_bundle.sbatch
and the queue item pins its digest, matching the dispatcher-only registration
route used by the PeakVI cross-cohort revision.  Nothing else about the bundle
differs from v9: same runner, same resources, same twenty-five run grid.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from masld_bench.artifacts import (
    freeze_tree,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan, source_lock


CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v10"
BUNDLE_ID = "cobolt-5seed-25run-v10-generic7"
SUPERSEDED_CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v9"
SUPERSEDED_BUNDLE_ID = "cobolt-5seed-25run-v9-generic6"
SUPERSEDES_FAILED_JOB_IDS = [21100294]
CONTINUATION_POLICY = {
    "startup": "verify_all_existing_run_roots_before_any_new_execution",
    "skip_only": (
        "attempt-001 verified succeeded by verify_run_execution_attempt "
        "and verify_frozen_tree"
    ),
    "failed_incomplete_or_ambiguous": "abort_before_any_new_execution",
    "mid_run_resume": False,
    "interrupted_run_recovery": "new_campaign_revision_required",
    "bundle_retry_identity": "distinct_immutable_bundle_and_queue_item_required",
}
SEEDS = (1103, 2107, 3109, 4111, 5113)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--superseded-candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--bundle-id", default=BUNDLE_ID)
    arguments = parser.parse_args()

    candidate = arguments.candidate.resolve(strict=True)
    superseded = arguments.superseded_candidate.resolve(strict=True)
    runner = arguments.runner.resolve(strict=True)
    if arguments.bundle_id != BUNDLE_ID:
        raise SystemExit("unexpected Cobolt v10 bundle revision")
    verify_frozen_tree(candidate)
    verify_frozen_tree(superseded)
    verify_frozen_tree(runner.parent)
    plan = load_frozen_plan(candidate)
    if plan["campaign"]["campaign_id"] != CAMPAIGN_ID:
        raise SystemExit("unexpected Cobolt campaign revision")

    superseded_plan = load_frozen_plan(superseded)
    if superseded_plan["campaign"]["campaign_id"] != SUPERSEDED_CAMPAIGN_ID:
        raise SystemExit("unexpected superseded Cobolt campaign revision")
    if plan["plan_sha256"] == superseded_plan["plan_sha256"]:
        raise SystemExit("v10 plan is not distinct from the superseded v9 plan")
    if plan["source_lock"] != superseded_plan["source_lock"]:
        raise SystemExit("v10 revision scope is not limited to the firewall baseline")
    if (
        plan["resource_firewall"]["snapshot_sha256"]
        == superseded_plan["resource_firewall"]["snapshot_sha256"]
    ):
        raise SystemExit("v10 did not re-baseline the Resource firewall snapshot")

    package_root = Path(str(plan["package_root"])).resolve(strict=True)
    included_paths = plan["source_lock"].get("included_paths")
    if source_lock(package_root, included_paths=included_paths) != plan["source_lock"]:
        raise SystemExit("source tree differs from the v10 candidate")

    runs = sorted(
        (
            run
            for run in plan["runs"]
            if run.get("model_id") == "cobolt"
            and run.get("resource_profile") == "gpu_rna_atac_smoke"
        ),
        key=lambda run: (int(run["seed"]), int(run["fold"])),
    )
    expected = [(seed, fold) for seed in SEEDS for fold in range(5)]
    if len(runs) != 25 or [
        (int(run["seed"]), int(run["fold"])) for run in runs
    ] != expected:
        raise SystemExit("Cobolt v10 must contain the frozen five-seed by five-fold grid")

    # The 25 v9 run roots must not be re-entered.  Run identities hash the
    # campaign name, so a genuine revision produces a disjoint set; assert it
    # rather than trusting it.
    superseded_ids = {str(run["run_id"]) for run in superseded_plan["runs"]}
    if superseded_ids & {str(run["run_id"]) for run in plan["runs"]}:
        raise SystemExit("v10 reuses superseded run identities")
    execution_output_root = candidate.parent / "executions" / candidate.name
    if execution_output_root.exists():
        raise SystemExit(f"v10 execution root already exists: {execution_output_root}")

    output = arguments.output
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": "masld-bench-sequential-run-bundle-v2",
        "bundle_id": arguments.bundle_id,
        "candidate": candidate.as_posix(),
        "candidate_manifest_sha256": sha256_file(candidate / "ARTIFACTS.json"),
        "plan_sha256": plan["plan_sha256"],
        "execution_output_root": execution_output_root.as_posix(),
        "bundle_receipt_root": (
            package_root / "executions" / "cobolt-5seed-bundles"
        ).as_posix(),
        "expected_resource_profile": "gpu_rna_atac_smoke",
        "expected_runtime_id": "gpu_rna_atac_torch_smoke",
        "expected_model_id": "cobolt",
        "expected_task_id": "rna_conditioned_atac",
        "expected_actions": ["prepare", "fit", "predict"],
        "continuation_policy": CONTINUATION_POLICY,
        "resources": {
            "account": "nslab",
            "partition": "gpu",
            "qos": "nslab",
            "gpu_type": "l40s",
            "gpus": 1,
            "cpus": 4,
            "memory_gb": 64,
            "walltime": "06:00:00",
        },
        "runner": runner.as_posix(),
        "runner_sha256": sha256_file(runner),
        "runs": [
            {
                "index": index,
                "seed": run["seed"],
                "fold": run["fold"],
                "run_id": run["run_id"],
            }
            for index, run in enumerate(runs, start=1)
        ],
    }
    write_json_exclusive(output / "bundle_manifest.json", manifest)
    write_json_exclusive(
        output / "supersedes.json",
        {
            "schema_version": "masld-bench-bundle-supersedes-v1",
            "bundle_id": arguments.bundle_id,
            "campaign_id": CAMPAIGN_ID,
            "superseded_bundle_id": SUPERSEDED_BUNDLE_ID,
            "superseded_campaign_id": SUPERSEDED_CAMPAIGN_ID,
            "supersedes_failed_job_ids": SUPERSEDES_FAILED_JOB_IDS,
            "superseded_candidate": superseded.as_posix(),
            "superseded_candidate_manifest_sha256": sha256_file(
                superseded / "ARTIFACTS.json"
            ),
            "superseded_plan_sha256": superseded_plan["plan_sha256"],
            "superseded_resource_snapshot_sha256": superseded_plan["resource_firewall"][
                "snapshot_sha256"
            ],
            "resource_snapshot_sha256": plan["resource_firewall"]["snapshot_sha256"],
            "revision_scope": (
                "Re-freeze the protected MASLD Resource firewall baseline against the "
                "current stable authority snapshot. Seeds, folds, model, adapter, "
                "hyperparameters, dataset view, donor split, runtime, resource profile, "
                "evaluator, and endpoint are unchanged."
            ),
            "mid_run_resume": False,
            "superseded_run_roots_reused": False,
            "superseded_run_roots_deleted": False,
            "incident": (
                "config/artifacts/incidents/"
                "cobolt_v9_resource_authority_firewall_21100294.json"
            ),
        },
    )
    process = subprocess.run(
        [
            sys.executable,
            runner.as_posix(),
            "--manifest",
            (output / "bundle_manifest.json").as_posix(),
            "--preflight",
        ],
        cwd=package_root,
        env={**os.environ, "PYTHONPATH": (package_root / "src").as_posix()},
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    write_text_exclusive(output / "preflight.stdout.txt", process.stdout)
    write_text_exclusive(output / "preflight.stderr.txt", process.stderr)
    if process.returncode != 0:
        raise SystemExit(f"bundle preflight failed: {process.stderr}")
    digest = freeze_tree(
        output,
        {
            "artifact_class": "gpu_sequential_dispatch_bundle",
            "bundle_id": arguments.bundle_id,
            "candidate_manifest_sha256": manifest["candidate_manifest_sha256"],
            "gpu_submission_authorized": False,
            "manual_gpu_sbatch_allowed": False,
            "walltime": "06:00:00",
            "continuation_policy": CONTINUATION_POLICY,
            "supersedes_failed_job_ids": SUPERSEDES_FAILED_JOB_IDS,
        },
    )
    print(
        json.dumps(
            {
                "bundle": output.resolve().as_posix(),
                "artifacts_sha256": digest,
                "bundle_manifest_sha256": sha256_file(output / "bundle_manifest.json"),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
