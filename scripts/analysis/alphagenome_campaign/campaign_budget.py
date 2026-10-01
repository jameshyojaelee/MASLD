#!/usr/bin/env python3
"""Account for individual campaign jobs and enforce remaining resource limits.

Submission uses an exclusive journal lock. Unfinished jobs retain their entire
reservation; terminal jobs release unused time only after SLURM accounting.
Existing P4 retrieval and unrelated allocations are excluded from this journal.
The user removed compute-hour and concurrent-GPU limits on 2026-09-16;
storage accounting and scheduler memory/partition requirements remain.
"""
from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[3]
DEFAULT = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/budget.json"
LIMITS = {"gpu_hours": None, "cpu_core_hours": None, "new_bytes": 1_500_000_000_000,
          "simultaneous_gpus": None}
TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}


def seconds(value):
    days, _, clock = value.rpartition("-")
    parts = list(map(int, clock.split(":")))
    if len(parts) == 1:
        parts = [0, parts[0], 0]  # Slurm's bare time is minutes
    elif len(parts) == 2:
        parts = [0, *parts]
    if len(parts) != 3 or min(parts) < 0:
        raise ValueError(f"Invalid wall time: {value}")
    return int(days or 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def parse_header(path):
    options = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#SBATCH "):
            continue
        fields = shlex.split(line[len("#SBATCH "):])
        if len(fields) == 1 and "=" in fields[0]:
            key, value = fields[0].split("=", 1)
        elif len(fields) == 2:
            key, value = fields
        else:
            raise ValueError(f"Use one explicit long option per header line: {line}")
        if key in options:
            raise ValueError(f"Repeated option: {key}")
        options[key] = value
    allowed = {"--job-name", "--partition", "--account", "--qos", "--gres", "--cpus-per-task",
               "--mem", "--time", "--output", "--error", "--dependency", "--chdir", "--nodes", "--ntasks"}
    if set(options) - allowed:
        raise ValueError(f"Unsupported allocation options: {set(options) - allowed}")
    if options.get("--account") != "nslab" or options.get("--qos") != "nslab":
        raise ValueError("This campaign requires account=nslab and qos=nslab")
    if options.get("--nodes", "1") != "1" or options.get("--ntasks", "1") != "1":
        raise ValueError("Only individual single-node, single-task jobs are supported")
    cpu = int(options["--cpus-per-task"])
    wall = seconds(options["--time"])
    if cpu <= 0 or wall <= 0:
        raise ValueError("Positive CPU and wall time required")
    gres = options.get("--gres", "")
    if gres and not re.fullmatch(r"gpu:[A-Za-z0-9_-]+:[1-9][0-9]*", gres):
        raise ValueError("GPU hardware and count must be explicit")
    gpu = int(gres.rsplit(":", 1)[-1]) if gres else 0
    partition = options.get("--partition")
    if partition not in {"cpu", "gpu", "bigmem", "io"} or (gpu > 0) != (partition == "gpu"):
        raise ValueError("Partition must match workload")
    if wall > {"cpu": 14, "gpu": 30, "bigmem": 7, "io": 3}[partition] * 86400:
        raise ValueError("Partition wall-time limit exceeded")
    memory = memory_bytes(options["--mem"])
    if memory <= 0 or (partition == "cpu" and (cpu > 16 or memory > 210 * 1024**3)):
        raise ValueError("Allocation exceeds CPU partition campaign rules")
    return options, {"cpus": cpu, "gpus": gpu, "wall_seconds": wall,
                     "memory_bytes": memory, "cpu_core_hours": cpu * wall / 3600,
                     "gpu_hours": gpu * wall / 3600}


def memory_bytes(value):
    match = re.fullmatch(r"([\d.]+)([KMGT]?)([cn]?)", value, re.I)
    if not match:
        raise ValueError(f"Unknown memory unit: {value}")
    return int(float(match[1]) * 1024 ** {"K": 1, "M": 2, "G": 3, "T": 4, "": 2}[match[2].upper()])


def reconcile(journal):
    ids = [str(j["job_id"]) for j in journal["jobs"] if j.get("job_id")]
    if not ids:
        return
    result = subprocess.run(["sacct", "-X", "-nP", "-j", ",".join(ids),
                             "--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,AllocTRES"],
                            check=True, capture_output=True, text=True)
    records = {r[0]: r for r in (line.split("|") for line in result.stdout.splitlines()) if len(r) >= 5}
    for job in journal["jobs"]:
        row = records.get(str(job.get("job_id")))
        if row is None:
            continue  # Unknown accounting is not zero usage.
        job["state"] = row[1].split()[0].rstrip("+")
        elapsed = int(row[2])
        job["elapsed_seconds"] = elapsed
        actual_cpu = int(row[3])
        gpu_match = re.search(r"(?:^|,)gres/gpu=(\d+)(?:,|$)", row[4])
        actual_gpu = int(gpu_match[1]) if gpu_match else job["gpus"]
        job["actual_cpu_core_hours"] = actual_cpu * elapsed / 3600
        job["actual_gpu_hours"] = actual_gpu * elapsed / 3600


def committed(journal, field):
    return sum(j.get("actual_" + field, j[field]) if j.get("state") in TERMINAL
               else j[field] for j in journal["jobs"])


# Archives this campaign created outside its own results root. The recorded
# 2026-09-17 amendment counts every deposit against one number, and these two
# retrievals are campaign deposits that happen to live under the Atlas tree.
# The pre-campaign P4 retrieval (p4-saturation-20260909T191917Z, jobs 21611396 /
# 21611739 / 21611740) stays the one exclusion, so the roots are named
# explicitly rather than globbed.
CAMPAIGN_ARCHIVES = (
    "GWAS/finemapping/results/alphagenome_atlas/p4-f5-background-20260917T145358Z",
    "GWAS/finemapping/results/alphagenome_atlas/p4-u2c-genetic-20260917T151313Z",
)


def _tree_bytes(root):
    if not root.exists():
        return 0
    return sum(p.stat().st_size for p in root.rglob("*") if p.is_file() and not p.is_symlink())


def archive_bytes(root, journal=None, refresh=False):
    """Sizes of the campaign archives under the Atlas tree, per path.

    These retrievals are complete and no longer written to, and walking a quarter
    of a million region files on every submission is a needless delay, so the
    measured sizes are kept in the journal and refreshed by `status`.
    """
    cached = (journal or {}).get("campaign_archive_bytes") or {}
    if not refresh and set(cached) == set(CAMPAIGN_ARCHIVES):
        return cached
    repo = root.resolve().parents[4]
    return {rel: _tree_bytes(repo / rel) for rel in CAMPAIGN_ARCHIVES}


def storage_bytes(root, journal=None, refresh=False):
    """Campaign results root plus the campaign archives deposited under the Atlas tree."""
    return _tree_bytes(root) + sum(archive_bytes(root, journal, refresh).values())


def apply_authorized_limits(journal):
    """Preserve the old envelope while recording the explicit user amendment."""
    if journal.get("limits") != LIMITS:
        journal.setdefault("resource_authorizations", []).append({
            "recorded_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "previous_limits": dict(journal.get("limits", {})),
            "limits": dict(LIMITS),
            "authority": "User instruction 2026-09-16: no more cap on compute; use needed CPU, GPU, RAM and cpu/bigmem/io/gpu partitions. User instruction 2026-09-17: raise the file allowance to 1.5 TB rather than re-scope what it counts, so every deposit keeps counting against one number.",
            "reservations_counted": "all deposits, including saturation archives under alphagenome_atlas/ and reference data under data/external/, even though the campaign results tree itself measured 10.3 GB when the allowance was raised; the pre-campaign P4 retrieval remains the one exclusion",
            "unchanged": "nslab account/QOS; individual jobs; scheduler and QOS memory limits; scientific holdouts, including the closed fold 0",
        })
        journal["limits"] = dict(LIMITS)


def check_compute_limits(journal, reservation):
    for field in ("gpu_hours", "cpu_core_hours"):
        limit = journal["limits"][field]
        if limit is not None and committed(journal, field) + reservation[field] > limit:
            raise ValueError(f"Insufficient remaining {field}")
    limit = journal["limits"]["simultaneous_gpus"]
    active_gpu = sum(j["gpus"] for j in journal["jobs"] if j.get("state") not in TERMINAL)
    if limit is not None and active_gpu + reservation["gpus"] > limit:
        raise ValueError("Concurrent GPU limit including pending reservations")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--journal", type=Path, default=DEFAULT)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status")
    submit = sub.add_parser("submit")
    submit.add_argument("--script", type=Path, required=True)
    submit.add_argument("--label", required=True)
    submit.add_argument("--storage-gb", type=float, required=True, help="Peak additional decimal GB")
    args = parser.parse_args()
    args.journal.parent.mkdir(parents=True, exist_ok=True)
    with args.journal.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        journal = json.loads(args.journal.read_text()) if args.journal.exists() else {
            "schema": "masld-week1-budget-v1", "limits": LIMITS, "jobs": [],
            "old_p4_retrieval_excluded": [21611396, 21611739, 21611740],
            "protected_outcomes": "closed",
        }
        apply_authorized_limits(journal)
        reconcile(journal)
        if args.command == "submit":
            if any(j["label"] == args.label for j in journal["jobs"]):
                raise ValueError("Label already recorded; review its disposition before submitting a new attempt")
            options, reservation = parse_header(args.script)
            check_compute_limits(journal, reservation)
            if args.storage_gb < 0:
                raise ValueError("Storage reservation must be nonnegative")
            reserved_bytes = sum(j["reserved_bytes"] for j in journal["jobs"])
            new_bytes = int(args.storage_gb * 1e9)
            if max(reserved_bytes, storage_bytes(args.journal.parent, journal)) + new_bytes > LIMITS["new_bytes"]:
                raise ValueError("Campaign storage budget exceeded")
            # squeue reports minimum memory per node (n) or per CPU (c).
            rows = subprocess.run(["squeue", "-h", "-u", __import__("getpass").getuser(),
                                   "-t", "RUNNING", "-o", "%m|%C|%D"],
                                  check=True, text=True, capture_output=True).stdout.splitlines()
            running_memory = 0
            for row in rows:
                mem, cpus, nodes = row.split("|")
                running_memory += memory_bytes(mem.strip()) * int(cpus if mem.strip().endswith("c") else nodes)
            if running_memory + reservation["memory_bytes"] > 7 * 1024**4:
                raise ValueError("User running-memory total plus request exceeds 7 TiB")
            print("\n".join(f"#SBATCH {k}={v}" for k, v in options.items()), flush=True)
            snapshot = args.journal.parent / "submitted_scripts" / (args.label + ".sbatch")
            snapshot.parent.mkdir(exist_ok=True)
            with snapshot.open("x") as handle:
                handle.write(args.script.read_text())
            entry = {"label": args.label, **reservation, "reserved_bytes": new_bytes,
                     "script": str(snapshot), "sha256": hashlib.sha256(snapshot.read_bytes()).hexdigest(),
                     "state": "SUBMITTING", "submitted_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
            journal["jobs"].append(entry)
            args.journal.write_text(json.dumps(journal, indent=2) + "\n")
            result = subprocess.run(["sbatch", "--parsable", str(snapshot)], capture_output=True, text=True)
            if result.returncode:
                entry.update(state="FAILED", actual_cpu_core_hours=0, actual_gpu_hours=0,
                             submission_error=result.stderr)
            else:
                entry.update(job_id=int(result.stdout.strip().split(";")[0]), state="PENDING")
            args.journal.write_text(json.dumps(journal, indent=2) + "\n")
            if result.returncode:
                raise RuntimeError(result.stderr)
        journal["campaign_archive_bytes"] = archive_bytes(args.journal.parent, journal, refresh=True)
        journal["measured_campaign_bytes"] = storage_bytes(args.journal.parent, journal)
        journal["accounted_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        args.journal.write_text(json.dumps(journal, indent=2) + "\n")
        print(json.dumps(journal, indent=2))


if __name__ == "__main__":
    main()
