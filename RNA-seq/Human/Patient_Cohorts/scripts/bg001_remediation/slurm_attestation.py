#!/usr/bin/env python3
"""Validate and capture authoritative SLURM execution evidence for a recount donor."""

from __future__ import annotations

import csv
import hashlib
import os
import re
import shlex
import shutil
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path


LEDGER_FIELDS = (
    "stage",
    "job_id",
    "dependency",
    "job_name",
    "partition",
    "qos",
    "array",
    "cpus_per_task",
    "memory",
    "time_limit",
    "submitted_utc",
)
SACCT_FIELDS = (
    "JobID",
    "JobIDRaw",
    "Cluster",
    "User",
    "UID",
    "Account",
    "AssocID",
    "DBIndex",
    "JobName",
    "State",
    "ExitCode",
    "Partition",
    "QOS",
    "ReqCPUS",
    "AllocCPUS",
    "ReqMem",
    "Timelimit",
    "Submit",
    "Start",
    "End",
    "WorkDir",
    "NodeList",
    "StdOut",
    "StdErr",
    "SubmitLine",
)
STAGE_SPECS = (
    {
        "stage": "bam_manifest_hash",
        "dependency_stage": "",
        "job_name": "samtools",
        "partition": "io",
        "qos": "interactive",
        "array": "0-7%4",
        "tasks": 8,
        "cpus": 1,
        "memory": "8G",
        "wrapper": "bam_manifest_hash_array.sbatch",
        "log_stem": "slurm_hash_%A_%a",
        "executed": True,
    },
    {
        "stage": "bam_manifest_finalize",
        "dependency_stage": "bam_manifest_hash",
        "job_name": "samtools",
        "partition": "io",
        "qos": "interactive",
        "array": "",
        "tasks": 1,
        "cpus": 1,
        "memory": "8G",
        "wrapper": "bam_manifest_finalize.sbatch",
        "log_stem": "slurm_finalize_%j",
        "executed": True,
    },
    {
        "stage": "freeze_baseline",
        "dependency_stage": "bam_manifest_finalize",
        "job_name": "edgeR",
        "partition": "io",
        "qos": "interactive",
        "array": "",
        "tasks": 1,
        "cpus": 2,
        "memory": "32G",
        "wrapper": "freeze_baseline.sbatch",
        "log_stem": "slurm_baseline_%j",
        "executed": True,
    },
    {
        "stage": "fragment_recount",
        "dependency_stage": "freeze_baseline",
        "job_name": "featureCounts",
        "partition": "io",
        "qos": "interactive",
        "array": "0-7%4",
        "tasks": 8,
        "cpus": 8,
        "memory": "32G",
        "wrapper": "recount_array.sbatch",
        "log_stem": "slurm_recount_%A_%a",
        "executed": True,
    },
    {
        "stage": "merge_and_validate",
        "dependency_stage": "fragment_recount",
        "job_name": "featureCounts",
        "partition": "io",
        "qos": "interactive",
        "array": "",
        "tasks": 1,
        "cpus": 4,
        "memory": "32G",
        "wrapper": "merge_validate.sbatch",
        "log_stem": "slurm_merge_%j",
        "executed": True,
    },
    {
        "stage": "four_arm_analysis",
        "dependency_stage": "merge_and_validate",
        "job_name": "limma",
        "partition": "cpu",
        "qos": "interactive",
        "array": "",
        "tasks": 1,
        "cpus": 8,
        "memory": "64G",
        "wrapper": "analysis_arms.sbatch",
        "log_stem": "slurm_analysis_%j",
        "executed": False,
    },
)
TIME_LIMIT = "2-00:00:00"
LEDGER_TIME_LIMIT = "48:00:00"
SUBMIT_TIME_TOLERANCE = timedelta(seconds=60)
EXPECTED_CLUSTER = "ne1"
EXPECTED_USER = "jameslee"
EXPECTED_ACCOUNT = "nslab"
EXPECTED_ASSOC_ID = "331"
RUN_RE = re.compile(r"^bg001-fragment-v211-gencode49-([0-9]{8})T[0-9]{6}Z$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def require_regular(path: Path, label: str) -> None:
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise SystemExit(f"Missing, empty, or symlinked {label}: {path}")


def parse_timestamp(value: str, label: str) -> datetime:
    if value in {"", "None", "Unknown"}:
        raise SystemExit(f"Missing scheduler timestamp for {label}")
    try:
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc)
    except ValueError as exc:
        raise SystemExit(f"Malformed scheduler timestamp for {label}: {value}") from exc


def read_submission_ledger(root: Path) -> list[dict[str, str]]:
    ledger = root / "contract/slurm_submission.tsv"
    require_regular(ledger, "source SLURM submission ledger")
    with ledger.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != LEDGER_FIELDS:
            raise SystemExit("Source SLURM submission ledger schema differs")
        rows = list(reader)
    if len(rows) != len(STAGE_SPECS):
        raise SystemExit("Source SLURM submission ledger must contain exactly six stages")
    by_stage = {row["stage"]: row for row in rows}
    if len(by_stage) != len(rows) or set(by_stage) != {spec["stage"] for spec in STAGE_SPECS}:
        raise SystemExit("Source SLURM submission ledger stages are missing or duplicated")
    ordered: list[dict[str, str]] = []
    job_ids: set[str] = set()
    for spec in STAGE_SPECS:
        row = by_stage[str(spec["stage"])]
        dependency_stage = str(spec["dependency_stage"])
        expected_dependency = (
            f"afterok:{by_stage[dependency_stage]['job_id']}" if dependency_stage else ""
        )
        expected = {
            "dependency": expected_dependency,
            "job_name": str(spec["job_name"]),
            "partition": str(spec["partition"]),
            "qos": str(spec["qos"]),
            "array": str(spec["array"]),
            "cpus_per_task": str(spec["cpus"]),
            "memory": str(spec["memory"]),
            "time_limit": LEDGER_TIME_LIMIT,
        }
        if any(row[field] != value for field, value in expected.items()):
            raise SystemExit(f"Source SLURM ledger resources/dependency differ for {spec['stage']}")
        if not row["job_id"].isdigit() or row["job_id"] in job_ids:
            raise SystemExit("Source SLURM ledger job IDs are malformed or duplicated")
        job_ids.add(row["job_id"])
        try:
            datetime.strptime(row["submitted_utc"], "%Y-%m-%dT%H:%M:%SZ")
        except ValueError as exc:
            raise SystemExit("Source SLURM ledger submitted_utc is malformed") from exc
        ordered.append(row)
    return ordered


def parse_sacct(raw: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for number, line in enumerate(raw.splitlines(), start=1):
        if not line:
            continue
        fields = next(csv.reader([line], delimiter="|"))
        if len(fields) != len(SACCT_FIELDS):
            raise SystemExit(f"Malformed sacct row {number}: expected {len(SACCT_FIELDS)} fields")
        rows.append(dict(zip(SACCT_FIELDS, fields)))
    return rows


def expected_submit_tokens(root: Path, row: dict[str, str], spec: dict[str, object]) -> list[str]:
    logs = root / "logs"
    remediation = root / "source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
    tokens = [
        "sbatch",
        "--parsable",
        "--no-requeue",
        f"--chdir={root}",
        f"--export=ALL,BG001_RUN_ROOT={root}",
    ]
    if row["dependency"]:
        tokens.append(f"--dependency={row['dependency']}")
    tokens.extend(
        (
            f"--output={logs}/{spec['log_stem']}.out",
            f"--error={logs}/{spec['log_stem']}.err",
            str(remediation / str(spec["wrapper"])),
        )
    )
    return tokens


def validate_empty_legacy_analysis(root: Path) -> None:
    for arm in ("R0", "F_locked", "F_legacy", "F_five"):
        arm_root = root / "arms" / arm
        if not arm_root.is_dir() or arm_root.is_symlink():
            raise SystemExit(f"Legacy source arm is missing or symlinked: {arm}")
        unexpected = [
            path
            for path in arm_root.rglob("*")
            if path.is_file() and path.name != ".bg001_candidate_root"
        ]
        if unexpected:
            raise SystemExit(f"Cancelled legacy analysis nevertheless wrote an arm artifact: {unexpected[0]}")
    if any((root / "logs").glob("slurm_analysis_*.out")) or any(
        (root / "logs").glob("slurm_analysis_*.err")
    ):
        raise SystemExit("Cancelled legacy analysis unexpectedly produced scheduler logs")


def observed_array_concurrency(
    normalized: list[dict[str, str]],
) -> dict[str, int]:
    result: dict[str, int] = {}
    for spec in STAGE_SPECS:
        if int(spec["tasks"]) <= 1:
            continue
        stage = str(spec["stage"])
        events: list[tuple[datetime, int]] = []
        for record in normalized:
            if record["stage"] != stage:
                continue
            events.append((parse_timestamp(record["Start"], f"{record['JobID']} start"), 1))
            events.append((parse_timestamp(record["End"], f"{record['JobID']} end"), -1))
        # Treat intervals as half-open: a completion at the exact instant of a
        # subsequent start frees its slot before the successor is counted.
        active = 0
        maximum = 0
        for _, delta in sorted(events, key=lambda event: (event[0], event[1])):
            active += delta
            if active < 0:
                raise SystemExit(f"Source SLURM array chronology is malformed: {stage}")
            maximum = max(maximum, active)
        if active != 0 or maximum > 4:
            raise SystemExit(f"Source SLURM array exceeded the four-task cap: {stage}")
        result[stage] = maximum
    return result


def validate_accounting(
    root: Path,
    ledger_rows: list[dict[str, str]],
    accounting_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    by_stage = {row["stage"]: row for row in ledger_rows}
    expected_total = sum(int(spec["tasks"]) for spec in STAGE_SPECS)
    if len(accounting_rows) != expected_total:
        raise SystemExit(
            f"Source SLURM attestation expected {expected_total} logical jobs, found {len(accounting_rows)}"
        )
    if any(not row["JobIDRaw"].isdigit() for row in accounting_rows):
        raise SystemExit("Source SLURM attestation contains malformed raw job IDs")
    if len({row["JobIDRaw"] for row in accounting_rows}) != len(accounting_rows):
        raise SystemExit("Source SLURM attestation contains duplicate raw job IDs")
    if any(not row["DBIndex"].isdigit() for row in accounting_rows):
        raise SystemExit("Source SLURM attestation contains malformed database indices")
    if len({row["DBIndex"] for row in accounting_rows}) != len(accounting_rows):
        raise SystemExit("Source SLURM attestation contains duplicate database indices")

    normalized: list[dict[str, str]] = []
    stage_times: dict[str, tuple[list[datetime], list[datetime]]] = {}
    for spec in STAGE_SPECS:
        ledger = by_stage[str(spec["stage"])]
        task_count = int(spec["tasks"])
        expected_ids = (
            {f"{ledger['job_id']}_{index}" for index in range(task_count)}
            if task_count > 1
            else {ledger["job_id"]}
        )
        observed = [row for row in accounting_rows if row["JobID"] in expected_ids]
        if {row["JobID"] for row in observed} != expected_ids or len(observed) != task_count:
            raise SystemExit(f"Source SLURM array/cardinality differs for {spec['stage']}")
        starts: list[datetime] = []
        ends: list[datetime] = []
        for record in observed:
            static_expected = {
                "Cluster": EXPECTED_CLUSTER,
                "User": EXPECTED_USER,
                "UID": str(root.stat().st_uid),
                "Account": EXPECTED_ACCOUNT,
                "AssocID": EXPECTED_ASSOC_ID,
                "JobName": str(spec["job_name"]),
                "Partition": str(spec["partition"]),
                "QOS": str(spec["qos"]),
                "ReqCPUS": str(spec["cpus"]),
                "ReqMem": str(spec["memory"]),
                "Timelimit": TIME_LIMIT,
                "WorkDir": str(root),
                "StdOut": str(root / "logs" / f"{spec['log_stem']}.out"),
                "StdErr": str(root / "logs" / f"{spec['log_stem']}.err"),
            }
            if any(record[field] != value for field, value in static_expected.items()):
                raise SystemExit(f"Source SLURM accounting resources/paths differ for {record['JobID']}")
            if shlex.split(record["SubmitLine"]) != expected_submit_tokens(root, ledger, spec):
                raise SystemExit(f"Source SLURM SubmitLine differs from frozen wrapper contract: {record['JobID']}")
            scheduler_submit = parse_timestamp(
                record["Submit"], f"{record['JobID']} submit"
            )
            ledger_submit = datetime.strptime(
                ledger["submitted_utc"], "%Y-%m-%dT%H:%M:%SZ"
            ).replace(tzinfo=timezone.utc)
            if abs(scheduler_submit - ledger_submit) > SUBMIT_TIME_TOLERANCE:
                raise SystemExit(
                    f"Source SLURM ledger/scheduler submit times differ for {record['JobID']}"
                )
            if bool(spec["executed"]):
                if (
                    record["State"] != "COMPLETED"
                    or record["ExitCode"] != "0:0"
                    or record["AllocCPUS"] != str(spec["cpus"])
                    or record["NodeList"] in {"", "None", "Unknown", "None assigned"}
                ):
                    raise SystemExit(f"Source SLURM executed job did not complete cleanly: {record['JobID']}")
                starts.append(parse_timestamp(record["Start"], f"{record['JobID']} start"))
                ends.append(parse_timestamp(record["End"], f"{record['JobID']} end"))
                if not scheduler_submit <= starts[-1] < ends[-1]:
                    raise SystemExit(
                        f"Source SLURM executed-job chronology differs: {record['JobID']}"
                    )
            else:
                if (
                    not record["State"].startswith("CANCELLED")
                    or record["ExitCode"] != "0:0"
                    or record["AllocCPUS"] != "0"
                    or record["Start"] not in {"", "None", "Unknown"}
                    or record["NodeList"] not in {"", "None", "Unknown", "None assigned"}
                ):
                    raise SystemExit("Legacy analysis was not cancelled before execution")
                cancelled_end = parse_timestamp(
                    record["End"], f"{record['JobID']} cancellation"
                )
                if scheduler_submit > cancelled_end:
                    raise SystemExit("Legacy analysis was cancelled before it was submitted")
            normalized.append({"stage": str(spec["stage"]), **record})
        if bool(spec["executed"]):
            stage_times[str(spec["stage"])] = (starts, ends)

    expected_job_ids = {
        row["JobID"]
        for row in normalized
    }
    if len(expected_job_ids) != expected_total or expected_job_ids != {
        row["JobID"] for row in accounting_rows
    }:
        raise SystemExit("Source SLURM attestation contains an extra or duplicate logical job")

    ordering = (
        ("bam_manifest_hash", "bam_manifest_finalize"),
        ("bam_manifest_finalize", "freeze_baseline"),
        ("freeze_baseline", "fragment_recount"),
        ("fragment_recount", "merge_and_validate"),
    )
    for predecessor, successor in ordering:
        if max(stage_times[predecessor][1]) > min(stage_times[successor][0]):
            raise SystemExit(f"Source SLURM temporal dependency order differs: {predecessor}->{successor}")
    validate_empty_legacy_analysis(root)
    stage_order = {str(spec["stage"]): index for index, spec in enumerate(STAGE_SPECS)}
    normalized.sort(key=lambda row: (stage_order[row["stage"]], row["JobID"]))
    return normalized


def capture(root_raw: Path) -> dict[str, object]:
    root = root_raw.resolve(strict=True)
    sentinel = root / ".bg001_candidate_root"
    if not sentinel.is_file() or sentinel.is_symlink():
        raise SystemExit("Source SLURM attestation requires a marked candidate root")
    run_match = RUN_RE.fullmatch(root.name)
    if run_match is None:
        raise SystemExit("Source SLURM attestation run ID is malformed")
    ledger_path = root / "contract/slurm_submission.tsv"
    require_regular(ledger_path, "source SLURM submission ledger")
    ledger_sha256 = sha256(ledger_path)
    ledger_rows = read_submission_ledger(root)
    if sha256(ledger_path) != ledger_sha256:
        raise SystemExit("Source SLURM submission ledger changed while it was read")
    job_ids = ",".join(row["job_id"] for row in ledger_rows)
    run_date = datetime.strptime(run_match.group(1), "%Y%m%d").date()
    start_date = (run_date - timedelta(days=1)).isoformat()
    sacct_executable_raw = shutil.which("sacct")
    if not sacct_executable_raw:
        raise SystemExit("sacct is unavailable for source execution attestation")
    sacct_executable = Path(sacct_executable_raw).resolve(strict=True)
    require_regular(sacct_executable, "sacct executable")
    sacct_status = sacct_executable.stat()
    if sacct_status.st_uid != 0 or sacct_status.st_mode & 0o022:
        raise SystemExit("sacct executable is not root-owned and non-writable by group/other")
    query = [
        str(sacct_executable),
        "-j",
        job_ids,
        "--array",
        "-X",
        "-P",
        "-n",
        "--local",
        "--duplicates",
        "--starttime",
        start_date,
        f"--format={','.join(SACCT_FIELDS)}",
    ]
    scheduler_environment = {**os.environ, "TZ": "UTC", "LC_ALL": "C"}
    completed = subprocess.run(
        query,
        check=True,
        text=True,
        capture_output=True,
        env=scheduler_environment,
    )
    if completed.stderr:
        raise SystemExit(f"sacct emitted stderr during attestation: {completed.stderr.strip()}")
    raw = completed.stdout
    accounting_rows = parse_sacct(raw)
    normalized = validate_accounting(root, ledger_rows, accounting_rows)
    array_concurrency = observed_array_concurrency(normalized)
    version = subprocess.run(
        [str(sacct_executable), "--version"],
        check=True,
        text=True,
        capture_output=True,
        env=scheduler_environment,
    )
    if version.stderr or not version.stdout.strip():
        raise SystemExit("Could not capture an exact sacct version")
    if sha256(ledger_path) != ledger_sha256:
        raise SystemExit("Source SLURM submission ledger changed during attestation")
    validation: dict[str, object] = {
        "schema_version": "bg001-source-slurm-attestation-v1",
        "status": "PASS",
        "source_run_id": root.name,
        "source_run_root": str(root),
        "ledger_relative_path": "contract/slurm_submission.tsv",
        "ledger_sha256": ledger_sha256,
        "raw_sacct_sha256": sha256_bytes(raw.encode()),
        "sacct_fields": list(SACCT_FIELDS),
        "sacct_query": query,
        "sacct_environment": {"TZ": "UTC", "LC_ALL": "C"},
        "sacct_executable": str(sacct_executable),
        "sacct_executable_sha256": sha256(sacct_executable),
        "sacct_version": version.stdout.strip(),
        "logical_jobs": len(normalized),
        "completed_jobs": sum(row["State"] == "COMPLETED" for row in normalized),
        "cancelled_unstarted_jobs": sum(row["State"].startswith("CANCELLED") for row in normalized),
        "observed_array_max_concurrency": array_concurrency,
        "normalized_records": normalized,
        "captured_utc": datetime.now(timezone.utc).isoformat(),
    }
    return {"raw_sacct_psv": raw, "validation": validation}


if __name__ == "__main__":
    raise SystemExit("slurm_attestation.py is a frozen import helper, not a standalone writer")
