#!/usr/bin/env python3
"""Build and freeze the five-seed Cobolt v9 GPU dispatch bundle."""

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


CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v9"
BUNDLE_ID = "cobolt-5seed-25run-v9-generic5"
WRAPPER_NAME = "model_work_212.sbatch"
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--bundle-id", default=BUNDLE_ID)
    arguments = parser.parse_args()

    candidate = arguments.candidate.resolve(strict=True)
    runner = arguments.runner.resolve(strict=True)
    if arguments.bundle_id not in {BUNDLE_ID, "cobolt-5seed-25run-v9-generic6"}:
        raise SystemExit("unexpected Cobolt v9 bundle revision")
    verify_frozen_tree(candidate)
    verify_frozen_tree(runner.parent)
    plan = load_frozen_plan(candidate)
    if plan["campaign"]["campaign_id"] != CAMPAIGN_ID:
        raise SystemExit("unexpected Cobolt campaign revision")
    package_root = Path(str(plan["package_root"])).resolve(strict=True)
    included_paths = plan["source_lock"].get("included_paths")
    if source_lock(package_root, included_paths=included_paths) != plan["source_lock"]:
        raise SystemExit("source tree differs from the v9 candidate")

    runs = sorted(
        (
            run
            for run in plan["runs"]
            if run.get("model_id") == "cobolt"
            and run.get("resource_profile") == "gpu_rna_atac_smoke"
        ),
        key=lambda run: (int(run["seed"]), int(run["fold"])),
    )
    expected = [
        (seed, fold)
        for seed in (1103, 2107, 3109, 4111, 5113)
        for fold in range(5)
    ]
    if len(runs) != 25 or [
        (int(run["seed"]), int(run["fold"])) for run in runs
    ] != expected:
        raise SystemExit("Cobolt v9 must contain the frozen five-seed by five-fold grid")

    output = arguments.output
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": "masld-bench-sequential-run-bundle-v2",
        "bundle_id": arguments.bundle_id,
        "candidate": candidate.as_posix(),
        "candidate_manifest_sha256": sha256_file(candidate / "ARTIFACTS.json"),
        "plan_sha256": plan["plan_sha256"],
        "execution_output_root": (
            candidate.parent / "executions" / candidate.name
        ).as_posix(),
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
    wrapper = f'''#!/usr/bin/env bash
#SBATCH --job-name=model-work-212
#SBATCH --account=nslab
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --gres=gpu:l40s:1
#SBATCH --qos=nslab
#SBATCH --output={package_root}/executions/model-work-212-%j.out
#SBATCH --error={package_root}/executions/model-work-212-%j.err

set -euo pipefail
umask 027
ROOT={package_root}
BUNDLE={output.resolve()}
CANDIDATE={candidate}
RUNNER={runner}
module load python/3.11.5-GCCcore-13.2.0
cd "${{ROOT}}"
PYTHONPATH=src python - "${{ROOT}}" "${{CANDIDATE}}" <<'PY'
from pathlib import Path
import sys
from masld_bench.planner import load_frozen_plan, source_lock
root = Path(sys.argv[1]).resolve(strict=True)
plan = load_frozen_plan(Path(sys.argv[2]))
included_paths = plan["source_lock"].get("included_paths")
if source_lock(root, included_paths=included_paths) != plan["source_lock"]:
    raise SystemExit("source tree differs from frozen v9 candidate; zero runs started")
PY
printf '%s  %s\\n' {sha256_file(runner)} "${{RUNNER}}" | sha256sum --check --strict
PYTHONPATH=src python "${{RUNNER}}" \\
  --manifest "${{BUNDLE}}/bundle_manifest.json" \\
  --execute --execution-id "${{SLURM_JOB_ID}}"
'''
    write_text_exclusive(output / WRAPPER_NAME, wrapper, mode=0o750)
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
            "walltime": "06:00:00",
            "continuation_policy": CONTINUATION_POLICY,
        },
    )
    print(
        json.dumps(
            {"bundle": output.resolve().as_posix(), "artifacts_sha256": digest},
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
