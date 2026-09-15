"""Readable view of the GPU rectangle: job -> bundle -> architecture.

Read-only. The real job names stay requirements-compliant with the dispatcher's
frozen `job_name_pattern`; this only renders them for a human.
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "executions/sequence-five-seed-rectangle-admission-21098556/plan/bundle_tasks.tsv"
LEDGER = ROOT / "executions/gpu-bundle-dispatch-state/submissions.tsv"
FRIENDLY = {
    "bpnet": "masld-bpnet",
    "chrombpnet": "masld-chrombpnet",
    "sequence_cnn_control": "masld-cnn-control",
    "sequence_transformer_control": "masld-transformer-control",
}


def main() -> int:
    bundle_model: dict[str, str] = {}
    planned: dict[str, int] = {}
    with PLAN.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            bundle_model.setdefault(row["bundle_id"], row["model_id"])
            planned[row["bundle_id"]] = planned.get(row["bundle_id"], 0) + 1

    job_bundle: dict[str, str] = {}
    for line in LEDGER.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) > 6:
            job_bundle[parts[6]] = parts[1]

    squeue = subprocess.run(
        ["squeue", "-h", "-o", "%i %t %M", "--me", "-p", "gpu"],
        capture_output=True, text=True, check=False,
    ).stdout.splitlines()

    print(f"{'job':<10}{'st':<4}{'elapsed':>10}  {'architecture':<28}{'bundle':<22}tasks")
    for line in sorted(squeue, key=lambda row: row.split()[1], reverse=True):
        job, state, elapsed = line.split()
        queue_bundle = job_bundle.get(job, "?")
        # Queue ids run +100 against the plan's bundle ids.
        plan_bundle = queue_bundle.replace("-8", "-7", 1) if queue_bundle.startswith("model-training-8") else queue_bundle
        model = bundle_model.get(plan_bundle, "?")
        total = planned.get(plan_bundle, 0)
        staged = list(ROOT.glob(f"executions/{plan_bundle}-{job}.staging/task_receipts/*.json"))
        done = len(staged) if staged else len(
            list(ROOT.glob(f"executions/{plan_bundle}-{job}/task_receipts/*.json"))
        )
        print(
            f"{job:<10}{state:<4}{elapsed:>10}  {FRIENDLY.get(model, model):<28}"
            f"{plan_bundle:<22}{done}/{total}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
