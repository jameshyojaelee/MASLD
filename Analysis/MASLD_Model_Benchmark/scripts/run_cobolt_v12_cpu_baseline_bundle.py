#!/usr/bin/env python3
"""Execute the frozen Cobolt v12 CPU control matrix in one bounded allocation.

Re-fits the 150 CPU baseline runs under campaign v12 because run_id is
canonical_sha256(identity_payload), which includes campaign_id and
code_lock_sha256; the B6K runtime revision re-keys every run in the campaign,
including these, whose own adapter did not change.  This supersedes
cobolt-v9-cpu-baselines--21088479 FOR v12 COMPARISON ONLY; the v9 matrix stays
intact as its own frozen record."""

from __future__ import annotations

import argparse
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile
import threading
import time
from typing import Any, Callable, Mapping, Sequence

from masld_bench.artifacts import (
    ArtifactRecord,
    freeze_tree,
    publish_directory_noreplace,
    verify_frozen_tree,
    write_json_exclusive,
)
from masld_bench.campaign import verify_run_execution_attempt
from masld_bench.contracts import ArtifactRef
from masld_bench.firewall import assert_resource_unchanged
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan, source_lock


CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v12"
PLAN_SHA256 = "3972f43112b8ecc4edf2aa5b128d0a8fcbc7953b21da9783ecd6ad302331c185"
CANDIDATE_SHA256 = "9b878eb08b3fc8a6012f473694b778a771da1bdd976f59e2966776f0bda3519f"
STABILITY_SHA256 = "0952f39b35c6d9ea437a32742bfacaf10d3a513e4507d6812b5ca6174e132e0a"
SOURCE_LOCK_SHA256 = "6b76abbe1264c3287682e8d395b1cd3677e59a551034fceab28c9f9d30e116df"
RESOURCE_SNAPSHOT_SHA256 = (
    "40c5f5c56e5b5695af2747b2c89b639b7d388c87f809aac712ae7d6e68c27c40"
)
CANDIDATE_BINDINGS = {
    CANDIDATE_SHA256: {
        "plan_sha256": PLAN_SHA256,
        "resource_snapshot_sha256": RESOURCE_SNAPSHOT_SHA256,
        "stability_sha256": STABILITY_SHA256,
        "stability_mode": "direct_candidate",
    },
}
EXPECTED_MODELS = (
    "assay_native_pseudobulk",
    "mean_track",
    "nearest_context",
    "shrunken_pseudobulk",
    "shuffled_context",
    "trans_only",
)
EXPECTED_SEEDS = (1103, 2107, 3109, 4111, 5113)
EXPECTED_FOLDS = tuple(range(5))
EXPECTED_RUNS = 150
EXPECTED_PROFILE = "cpu_baseline_smoke"
EXPECTED_TASK = "rna_conditioned_atac"
EXPECTED_DATASET = "gse296875"
EXPECTED_ACTIONS = ["prepare", "fit", "predict"]
RUN_CPUS = 1
RUN_MEMORY_GB = 16
MAX_WORKERS = 12
SHA256 = re.compile(r"^[0-9a-f]{64}$")
NESTED_SBATCH = re.compile(r"(^|[;&|]\s*)sbatch(?:\s|$)")
CONTINUATION_POLICY = {
    "startup": "verify_all_existing_run_roots_before_any_new_execution",
    "skip_only": (
        "attempt-001 verified succeeded by verify_run_execution_attempt "
        "and verify_frozen_tree"
    ),
    "failed_incomplete_or_ambiguous": "abort_before_any_new_execution",
    "mid_run_resume": False,
    "interrupted_run_recovery": "new_campaign_revision_required",
    "stop_launching": "first_new_run_failure_or_usr1",
    "bundle_retry_identity": "distinct_immutable_slurm_job_receipt_required",
}


class CpuBundleError(RuntimeError):
    """Raised when bundle execution would weaken a frozen requirement."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CpuBundleError(f"cannot read {label}: {error}") from error
    if not isinstance(value, dict):
        raise CpuBundleError(f"{label} must contain a JSON object")
    return value


def _assert_sha256(value: str, label: str) -> None:
    if not SHA256.fullmatch(value):
        raise CpuBundleError(f"{label} must be a lowercase SHA-256 digest")


def _verify_bound_tree(path: Path, expected_sha256: str, label: str) -> Path:
    _assert_sha256(expected_sha256, f"{label} binding")
    if not path.is_absolute():
        raise CpuBundleError(f"{label} path must be absolute")
    try:
        resolved = path.resolve(strict=True)
    except OSError as error:
        raise CpuBundleError(f"{label} path is unavailable: {path}") from error
    if path.is_symlink() or not resolved.is_dir():
        raise CpuBundleError(f"{label} must be a non-symlink directory")
    verify_frozen_tree(resolved)
    if sha256_file(resolved / "ARTIFACTS.json") != expected_sha256:
        raise CpuBundleError(f"{label} ARTIFACTS SHA-256 differs")
    return resolved


def _check_capacity(
    *, workers: int, allocated_cpus: int, allocated_memory_gb: int
) -> None:
    if workers < 1 or workers > MAX_WORKERS:
        raise CpuBundleError(f"workers must be between 1 and {MAX_WORKERS}")
    if allocated_cpus < workers * RUN_CPUS:
        raise CpuBundleError("allocation has fewer CPUs than the worker contract")
    if allocated_memory_gb < workers * RUN_MEMORY_GB:
        raise CpuBundleError("allocation has less memory than the worker contract")


def _validate_job_script(candidate: Path, run_id: str) -> str:
    script = candidate / "jobs" / f"{run_id}.sbatch"
    if not script.is_file() or script.is_symlink():
        raise CpuBundleError(f"frozen run script is unavailable: {script}")
    text = script.read_text(encoding="utf-8")
    body = "\n".join(line for line in text.splitlines() if not line.startswith("#"))
    if "#SBATCH --array" in text or NESTED_SBATCH.search(body):
        raise CpuBundleError(f"array or nested sbatch is forbidden: {script}")
    required = (
        "#SBATCH --partition=cpu",
        "#SBATCH --cpus-per-task=1",
        "#SBATCH --mem=16G",
        "#SBATCH --time=01:00:00",
        "-m masld_bench.cli run execute",
        f"--run-id {run_id}",
    )
    if any(value not in text for value in required):
        raise CpuBundleError(f"frozen run script contract differs: {script}")
    return sha256_file(script)


def select_cpu_runs(plan: Mapping[str, Any], candidate: Path) -> list[dict[str, Any]]:
    """Return the exact six-family, five-seed, five-fold CPU roster."""

    candidate_manifest_sha256 = sha256_file(candidate / "ARTIFACTS.json")
    binding = CANDIDATE_BINDINGS.get(candidate_manifest_sha256)
    if binding is None or plan.get("plan_sha256") != binding["plan_sha256"]:
        raise CpuBundleError("candidate plan SHA-256 differs from the reviewed v12 plan")
    campaign = plan.get("campaign")
    if not isinstance(campaign, Mapping) or campaign.get("campaign_id") != CAMPAIGN_ID:
        raise CpuBundleError("candidate campaign identity differs")
    if plan.get("source_lock_sha256") != SOURCE_LOCK_SHA256:
        raise CpuBundleError("candidate source-lock identity differs")
    firewall = plan.get("resource_firewall")
    if (
        not isinstance(firewall, Mapping)
        or firewall.get("snapshot_sha256") != binding["resource_snapshot_sha256"]
    ):
        raise CpuBundleError("candidate Resource snapshot identity differs")
    if plan.get("safety") != {
        "arrays": False,
        "sealed_features": False,
        "sealed_labels": False,
        "submit_enabled": True,
    }:
        raise CpuBundleError("candidate safety contract differs")
    if plan.get("retry_policy") != {
        "max_attempts": 1,
        "retryable_states": [],
        "checkpoint_resume_required": True,
        "software_correction": "new_campaign_revision",
        "reason": (
            "No standardized immutable resume-checkpoint ArtifactRef contract "
            "is implemented; repeated execution fails closed."
        ),
    }:
        raise CpuBundleError("candidate retry policy differs")

    raw_runs = plan.get("runs")
    if not isinstance(raw_runs, list):
        raise CpuBundleError("candidate has no run inventory")
    selected = sorted(
        (
            dict(run)
            for run in raw_runs
            if isinstance(run, Mapping)
            and run.get("resource_profile") == EXPECTED_PROFILE
        ),
        key=lambda run: str(run.get("run_id", "")),
    )
    if len(selected) != EXPECTED_RUNS:
        raise CpuBundleError(
            f"expected {EXPECTED_RUNS} CPU baseline runs, observed {len(selected)}"
        )
    model_counts = Counter(str(run.get("model_id", "")) for run in selected)
    if model_counts != Counter({model_id: 25 for model_id in EXPECTED_MODELS}):
        raise CpuBundleError("CPU baseline model counts differ")

    expected_grid = {
        (model_id, seed, fold)
        for model_id in EXPECTED_MODELS
        for seed in EXPECTED_SEEDS
        for fold in EXPECTED_FOLDS
    }
    observed_grid: set[tuple[str, int, int]] = set()
    for run in selected:
        run_id = str(run.get("run_id", ""))
        if not SHA256.fullmatch(run_id):
            raise CpuBundleError(f"invalid CPU run ID: {run_id}")
        try:
            key = (str(run["model_id"]), int(run["seed"]), int(run["fold"]))
        except (KeyError, TypeError, ValueError) as error:
            raise CpuBundleError(f"invalid CPU run grid entry: {run_id}") from error
        if key in observed_grid:
            raise CpuBundleError(f"duplicate CPU run grid entry: {key}")
        observed_grid.add(key)
        metadata = run.get("metadata")
        if not isinstance(metadata, Mapping):
            raise CpuBundleError(f"CPU run metadata is invalid: {run_id}")
        routing = metadata.get("dataset_routing")
        expected_routing = routing.get(EXPECTED_DATASET) if isinstance(routing, Mapping) else None
        if (
            run.get("stage") != "smoke"
            or run.get("task_id") != EXPECTED_TASK
            or run.get("dataset_ids") != [EXPECTED_DATASET]
            or run.get("action") != EXPECTED_ACTIONS
            or run.get("checkpoint") is not None
            or metadata.get("fit_dataset_ids") != [EXPECTED_DATASET]
            or metadata.get("development_prediction_dataset_ids") != []
            or metadata.get("prediction_first_stress_dataset_ids") != []
            or metadata.get("sealed_prediction_dataset_ids") != []
            or not isinstance(expected_routing, Mapping)
            or expected_routing.get("fit_allowed") is not True
            or expected_routing.get("prediction_first_outcome_artifact_included")
            is not False
        ):
            raise CpuBundleError(f"CPU run data firewall differs: {run_id}")
        run["job_script_sha256"] = _validate_job_script(candidate, run_id)
    if observed_grid != expected_grid:
        raise CpuBundleError("CPU run seed/fold grid differs")
    return selected


def _verify_unique_inputs(runs: Sequence[Mapping[str, Any]]) -> int:
    unique: dict[tuple[str, str, int, str | None], ArtifactRef] = {}
    for run in runs:
        raw_inputs = run.get("inputs")
        if not isinstance(raw_inputs, list):
            raise CpuBundleError(f"run has no input inventory: {run.get('run_id')}")
        for raw in raw_inputs:
            if not isinstance(raw, Mapping):
                raise CpuBundleError("run input inventory contains a non-object")
            artifact = ArtifactRef.from_dict(raw)
            key = (artifact.path, artifact.sha256, artifact.size_bytes, artifact.role)
            unique[key] = artifact
    for artifact in unique.values():
        artifact.validate()
    return len(unique)


def validate_startup(
    *,
    candidate_path: Path,
    candidate_sha256: str,
    stability_path: Path,
    stability_sha256: str,
    workers: int,
    allocated_cpus: int,
    allocated_memory_gb: int,
) -> tuple[Path, Path, dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    """Validate every shared input before any run root can be created."""

    _check_capacity(
        workers=workers,
        allocated_cpus=allocated_cpus,
        allocated_memory_gb=allocated_memory_gb,
    )
    binding = CANDIDATE_BINDINGS.get(candidate_sha256)
    if binding is None:
        raise CpuBundleError("candidate SHA-256 is not the reviewed v12 candidate")
    if stability_sha256 != binding["stability_sha256"]:
        raise CpuBundleError("stability SHA-256 is not the reviewed v12 audit")
    candidate = _verify_bound_tree(candidate_path, candidate_sha256, "candidate")
    plan = load_frozen_plan(candidate)
    runs = select_cpu_runs(plan, candidate)

    package_root = Path(str(plan.get("package_root", ""))).resolve(strict=True)
    included_paths = plan.get("source_lock", {}).get("included_paths")
    if source_lock(package_root, included_paths=included_paths) != plan.get("source_lock"):
        raise CpuBundleError("source tree differs from the reviewed v12 candidate")

    stability = _verify_bound_tree(stability_path, stability_sha256, "stability")
    stability_receipt = _read_json_object(
        stability / "receipt.json", "Resource stability receipt"
    )
    # Bound to the plan and candidate digests rather than to campaign_id: a digest
    # pins the exact reviewed output file, a name does not.
    if (
        stability_receipt.get("schema_version")
        != "masld-bench-cobolt-v12-revision-scope-validation-v1"
        or stability_receipt.get("plan_sha256") != PLAN_SHA256
        or stability_receipt.get("candidate_manifest_sha256") != CANDIDATE_SHA256
        or stability_receipt.get("resource_firewall_snapshot_at_freeze")
        != binding["resource_snapshot_sha256"]
        or stability_receipt.get("resource_firewall_still_current") is not True
        or stability_receipt.get("scientific_content_identical_to_v11") is not True
        or stability_receipt.get("l40s_validation_reused") is not False
        or stability_receipt.get("status") != "passed_firewall_baseline_only_revision"
    ):
        raise CpuBundleError("Resource stability receipt differs")

    repository_root = Path(str(plan.get("repository_root", ""))).resolve(strict=True)
    observed_firewall = assert_resource_unchanged(
        plan["resource_firewall"], repository_root=repository_root
    )
    if observed_firewall.get("snapshot_sha256") != binding["resource_snapshot_sha256"]:
        raise CpuBundleError("Resource authorities changed after stability validation")
    unique_inputs = _verify_unique_inputs(runs)
    startup = {
        "schema_version": "masld-bench-cobolt-v12-cpu-startup-v1",
        "campaign_id": CAMPAIGN_ID,
        "candidate": candidate.as_posix(),
        "candidate_manifest_sha256": candidate_sha256,
        "plan_sha256": binding["plan_sha256"],
        "source_lock_sha256": SOURCE_LOCK_SHA256,
        "stability": stability.as_posix(),
        "stability_manifest_sha256": stability_sha256,
        "resource_snapshot_sha256": binding["resource_snapshot_sha256"],
        "run_count": len(runs),
        "unique_input_artifacts_verified": unique_inputs,
        "workers": workers,
        "allocated_cpus": allocated_cpus,
        "allocated_memory_gb": allocated_memory_gb,
        "sealed_features_opened": False,
        "sealed_labels_opened": False,
        "outcome_records_parsed_or_joined": False,
        "status": "pass",
    }
    return candidate, stability, plan, runs, startup


def _attempt_record(
    attempt: Path,
    *,
    run_id: str,
    require_succeeded: bool,
    attempt_verifier: Callable[..., Mapping[str, Any]] = verify_run_execution_attempt,
    tree_verifier: Callable[[Path], Mapping[str, Any]] = verify_frozen_tree,
) -> dict[str, Any]:
    tree_verifier(attempt)
    receipt = attempt_verifier(attempt, require_succeeded=require_succeeded)
    if str(receipt.get("run_id", "")) != run_id:
        raise CpuBundleError(f"attempt receipt run identity differs: {run_id}")
    return {
        "attempt_path": attempt.as_posix(),
        "run_status": str(receipt.get("status", "unknown")),
        "run_receipt_sha256": sha256_file(attempt / "run_execution_receipt.json"),
        "run_artifacts_sha256": sha256_file(attempt / "ARTIFACTS.json"),
        "run_complete_sha256": sha256_file(attempt / "COMPLETE"),
    }


def classify_existing_successes(
    runs: Sequence[Mapping[str, Any]],
    execution_root: Path,
    *,
    campaign_run_ids: frozenset[str] | None = None,
    attempt_verifier: Callable[..., Mapping[str, Any]] = verify_run_execution_attempt,
    tree_verifier: Callable[[Path], Mapping[str, Any]] = verify_frozen_tree,
) -> dict[str, dict[str, Any]]:
    """Fail closed on any prior state other than one verified success."""

    verified: dict[str, dict[str, Any]] = {}
    errors: list[str] = []
    selected_ids = {str(run["run_id"]) for run in runs}
    # The foreign-run guard allows any run id belonging to THIS frozen campaign,
    # not only the subset this bundle executes. The evaluator resolves all 175
    # runs from a single execution root, so the GPU arm's 25 run roots must be
    # allowed to coexist here with the 150 CPU controls. Anything outside the
    # frozen plan is still rejected, which is what the guard exists for.
    permitted_ids = selected_ids if campaign_run_ids is None else (selected_ids | campaign_run_ids)
    runs_root = execution_root / "runs"
    if runs_root.exists():
        if runs_root.is_symlink() or not runs_root.is_dir():
            raise CpuBundleError("execution runs root is not a regular directory")
        unknown = sorted(
            child.name
            for child in runs_root.iterdir()
            if child.name not in permitted_ids
        )
        if unknown:
            raise CpuBundleError(
                "execution root contains unknown run identities: " + ", ".join(unknown[:3])
            )
    for run in runs:
        run_id = str(run["run_id"])
        run_root = runs_root / run_id
        if not run_root.exists():
            if run_root.is_symlink():
                errors.append(f"{run_id}: dangling run-root symlink")
            continue
        if run_root.is_symlink() or not run_root.is_dir():
            errors.append(f"{run_id}: run root is not a regular directory")
            continue
        children = sorted(run_root.iterdir(), key=lambda item: item.name)
        if [child.name for child in children] != ["attempt-001"]:
            errors.append(f"{run_id}: run root is incomplete or ambiguous")
            continue
        attempt = children[0]
        if attempt.is_symlink() or not attempt.is_dir():
            errors.append(f"{run_id}: attempt-001 is not a regular directory")
            continue
        try:
            record = _attempt_record(
                attempt,
                run_id=run_id,
                require_succeeded=True,
                attempt_verifier=attempt_verifier,
                tree_verifier=tree_verifier,
            )
        except Exception as error:
            errors.append(f"{run_id}: attempt-001 verification failed: {error}")
            continue
        if record["run_status"] != "succeeded":
            errors.append(f"{run_id}: attempt-001 is not successful")
            continue
        verified[run_id] = record
    if errors:
        raise CpuBundleError(
            "existing run roots failed startup verification; zero new runs started: "
            + "; ".join(errors[:3])
        )
    return verified


def _direct_command(
    *, plan: Mapping[str, Any], candidate: Path, execution_root: Path, run_id: str
) -> list[str]:
    control = plan.get("control_plane_lock")
    if not isinstance(control, Mapping):
        raise CpuBundleError("candidate control-plane lock is invalid")
    python = Path(str(control.get("python_executable", ""))).resolve(strict=True)
    if not python.is_file() or python.is_symlink():
        raise CpuBundleError("frozen control-plane Python is unavailable")
    return [
        python.as_posix(),
        "-m",
        "masld_bench.cli",
        "run",
        "execute",
        "--candidate",
        candidate.as_posix(),
        "--run-id",
        run_id,
        "--output-root",
        execution_root.as_posix(),
    ]


def _run_one(
    *,
    run: Mapping[str, Any],
    plan: Mapping[str, Any],
    candidate: Path,
    execution_root: Path,
    staging: Path,
    command_runner: Callable[..., subprocess.CompletedProcess[bytes]] = subprocess.run,
) -> dict[str, Any]:
    run_id = str(run["run_id"])
    command = _direct_command(
        plan=plan,
        candidate=candidate,
        execution_root=execution_root,
        run_id=run_id,
    )
    if any(part == "sbatch" for part in command):
        raise CpuBundleError("nested sbatch command construction is forbidden")
    stdout_path = staging / "logs" / f"{run_id}.out"
    stderr_path = staging / "logs" / f"{run_id}.err"
    started_at = _utc_now()
    monotonic_start = time.monotonic()
    control = plan["control_plane_lock"]
    environment = {
        **os.environ,
        "PYTHONPATH": str(control["pythonpath"]),
        "PYTHONNOUSERSITE": "1",
        "PYTHONDONTWRITEBYTECODE": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    }
    process_error: str | None = None
    returncode: int | None = None
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        try:
            process = command_runner(
                command,
                cwd=str(plan["package_root"]),
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
            returncode = process.returncode
        except Exception as error:
            process_error = f"{type(error).__name__}: {error}"

    attempt = execution_root / "runs" / run_id / "attempt-001"
    terminal: dict[str, Any] = {
        "attempt_path": attempt.as_posix() if attempt.exists() else None,
        "run_status": "missing_terminal_attempt",
        "run_receipt_sha256": None,
        "run_artifacts_sha256": None,
        "run_complete_sha256": None,
        "verification_error": None,
    }
    if attempt.is_dir() and not attempt.is_symlink():
        try:
            terminal.update(
                _attempt_record(
                    attempt,
                    run_id=run_id,
                    require_succeeded=False,
                )
            )
        except Exception as error:
            terminal["verification_error"] = f"{type(error).__name__}: {error}"
    else:
        terminal["verification_error"] = "attempt-001 is unavailable"
    disposition = (
        "succeeded"
        if returncode == 0
        and terminal["run_status"] == "succeeded"
        and terminal["verification_error"] is None
        else "failed"
    )
    return {
        "schema_version": "masld-bench-cobolt-v12-cpu-run-disposition-v1",
        "run_id": run_id,
        "model_id": str(run["model_id"]),
        "seed": int(run["seed"]),
        "fold": int(run["fold"]),
        "job_script_sha256": str(run["job_script_sha256"]),
        "job_script_invoked": False,
        "execution_mode": "direct_frozen_control_plane_no_nested_sbatch",
        "command": command,
        "started_at": started_at,
        "finished_at": _utc_now(),
        "elapsed_seconds": time.monotonic() - monotonic_start,
        "returncode": returncode,
        "process_error": process_error,
        "stdout": ArtifactRecord.from_path(stdout_path, relative_to=staging).as_dict(),
        "stderr": ArtifactRecord.from_path(stderr_path, relative_to=staging).as_dict(),
        **terminal,
        "disposition": disposition,
        "sealed_features_opened": False,
        "sealed_labels_opened": False,
        "outcome_records_parsed_or_joined": False,
    }


def _skipped_record(run: Mapping[str, Any], terminal: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "masld-bench-cobolt-v12-cpu-run-disposition-v1",
        "run_id": str(run["run_id"]),
        "model_id": str(run["model_id"]),
        "seed": int(run["seed"]),
        "fold": int(run["fold"]),
        "job_script_sha256": str(run["job_script_sha256"]),
        "job_script_invoked": False,
        "execution_mode": "reuse_verified_attempt_no_command",
        "command": None,
        "started_at": None,
        "finished_at": _utc_now(),
        "elapsed_seconds": 0.0,
        "returncode": None,
        "process_error": None,
        "stdout": None,
        "stderr": None,
        **dict(terminal),
        "verification_error": None,
        "disposition": "skipped_verified_success",
        "sealed_features_opened": False,
        "sealed_labels_opened": False,
        "outcome_records_parsed_or_joined": False,
    }


def _not_started_record(run: Mapping[str, Any], reason: str) -> dict[str, Any]:
    return {
        "schema_version": "masld-bench-cobolt-v12-cpu-run-disposition-v1",
        "run_id": str(run["run_id"]),
        "model_id": str(run["model_id"]),
        "seed": int(run["seed"]),
        "fold": int(run["fold"]),
        "job_script_sha256": str(run["job_script_sha256"]),
        "job_script_invoked": False,
        "execution_mode": "not_started",
        "command": None,
        "started_at": None,
        "finished_at": _utc_now(),
        "elapsed_seconds": 0.0,
        "returncode": None,
        "process_error": None,
        "stdout": None,
        "stderr": None,
        "attempt_path": None,
        "run_status": "not_started",
        "run_receipt_sha256": None,
        "run_artifacts_sha256": None,
        "run_complete_sha256": None,
        "verification_error": reason,
        "disposition": "not_started",
        "sealed_features_opened": False,
        "sealed_labels_opened": False,
        "outcome_records_parsed_or_joined": False,
    }


def _write_record(staging: Path, index: int, record: Mapping[str, Any]) -> None:
    write_json_exclusive(
        staging / "records" / f"{index:03d}-{record['run_id']}.json", record
    )


def _acquire_execution_lock(execution_root: Path) -> tuple[int, Path]:
    execution_root.mkdir(parents=True, exist_ok=True)
    if execution_root.is_symlink():
        raise CpuBundleError("execution root may not be a symlink")
    lock_path = execution_root / ".cpu-bundle.lock"
    flags = os.O_CREAT | os.O_RDWR
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(lock_path, flags, 0o640)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        os.close(descriptor)
        raise CpuBundleError("another CPU bundle owns the execution root") from error
    return descriptor, lock_path


def execute_bundle(arguments: argparse.Namespace) -> int:
    candidate, stability, plan, runs, startup = validate_startup(
        candidate_path=arguments.candidate,
        candidate_sha256=arguments.candidate_sha256,
        stability_path=arguments.stability,
        stability_sha256=arguments.stability_sha256,
        workers=arguments.workers,
        allocated_cpus=arguments.allocated_cpus,
        allocated_memory_gb=arguments.allocated_memory_gb,
    )
    if arguments.executor_sha256 != sha256_file(Path(__file__).resolve()):
        raise CpuBundleError("executor SHA-256 differs from the reviewed binding")
    if not re.fullmatch(r"[0-9]+", arguments.execution_id):
        raise CpuBundleError("execution ID must be the numeric SLURM job ID")
    execution_root = arguments.execution_root.resolve()
    receipt_root = arguments.receipt_root.resolve()
    if (
        execution_root == receipt_root
        or execution_root in receipt_root.parents
        or receipt_root in execution_root.parents
    ):
        raise CpuBundleError("execution and receipt roots may not contain one another")
    receipt_root.mkdir(parents=True, exist_ok=True)
    target = receipt_root / f"cobolt-v12-cpu-baselines--{arguments.execution_id}"
    if target.exists() or target.is_symlink():
        raise CpuBundleError(f"bundle receipt target already exists: {target}")

    lock_descriptor, lock_path = _acquire_execution_lock(execution_root)
    try:
        verified = classify_existing_successes(
            runs,
            execution_root,
            campaign_run_ids=frozenset(str(r["run_id"]) for r in plan["runs"]),
        )
        staging = Path(
            tempfile.mkdtemp(prefix=f".{target.name}.", dir=receipt_root)
        )
        (staging / "logs").mkdir()
        (staging / "records").mkdir()
        write_json_exclusive(staging / "startup_receipt.json", startup)
        stopped = threading.Event()

        def request_stop(signum: int, _frame: object) -> None:
            if signum == signal.SIGUSR1:
                stopped.set()

        previous_usr1 = signal.signal(signal.SIGUSR1, request_stop)
        records_by_id: dict[str, dict[str, Any]] = {}
        run_by_id = {str(run["run_id"]): run for run in runs}
        index_by_id = {
            str(run["run_id"]): index for index, run in enumerate(runs, start=1)
        }
        started_at = _utc_now()
        try:
            pending: deque[dict[str, Any]] = deque()
            for run in runs:
                run_id = str(run["run_id"])
                if run_id in verified:
                    record = _skipped_record(run, verified[run_id])
                    records_by_id[run_id] = record
                    _write_record(staging, index_by_id[run_id], record)
                else:
                    pending.append(run)

            with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
                active: dict[Future[dict[str, Any]], str] = {}

                def launch_available() -> None:
                    while pending and len(active) < arguments.workers and not stopped.is_set():
                        run = pending.popleft()
                        run_id = str(run["run_id"])
                        future = executor.submit(
                            _run_one,
                            run=run,
                            plan=plan,
                            candidate=candidate,
                            execution_root=execution_root,
                            staging=staging,
                        )
                        active[future] = run_id

                launch_available()
                while active:
                    completed, _ = wait(active, return_when=FIRST_COMPLETED)
                    for future in completed:
                        run_id = active.pop(future)
                        try:
                            record = future.result()
                        except Exception as error:
                            record = _not_started_record(
                                run_by_id[run_id],
                                f"worker_exception={type(error).__name__}: {error}",
                            )
                            record["disposition"] = "failed"
                            record["run_status"] = "worker_exception"
                        records_by_id[run_id] = record
                        _write_record(staging, index_by_id[run_id], record)
                        if record["disposition"] != "succeeded":
                            stopped.set()
                    launch_available()

            stop_reason = (
                "new_run_failure_or_usr1_stopped_further_launches"
                if stopped.is_set()
                else "not_applicable"
            )
            while pending:
                run = pending.popleft()
                run_id = str(run["run_id"])
                record = _not_started_record(run, stop_reason)
                records_by_id[run_id] = record
                _write_record(staging, index_by_id[run_id], record)
        finally:
            signal.signal(signal.SIGUSR1, previous_usr1)

        ordered = [records_by_id[str(run["run_id"])] for run in runs]
        counts = Counter(record["disposition"] for record in ordered)
        succeeded = counts["succeeded"] + counts["skipped_verified_success"]
        status = "succeeded" if succeeded == EXPECTED_RUNS else "partial_failure"
        receipt = {
            "schema_version": "masld-bench-cobolt-v12-cpu-bundle-receipt-v1",
            "campaign_id": CAMPAIGN_ID,
            "candidate": candidate.as_posix(),
            "candidate_manifest_sha256": arguments.candidate_sha256,
            "plan_sha256": plan["plan_sha256"],
            "source_lock_sha256": SOURCE_LOCK_SHA256,
            "stability": stability.as_posix(),
            "stability_manifest_sha256": arguments.stability_sha256,
            "resource_snapshot_sha256": plan["resource_firewall"]["snapshot_sha256"],
            "executor": Path(__file__).resolve().as_posix(),
            "executor_sha256": arguments.executor_sha256,
            "execution_root": execution_root.as_posix(),
            "execution_lock": lock_path.as_posix(),
            "execution_id": arguments.execution_id,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "started_at": started_at,
            "finished_at": _utc_now(),
            "status": status,
            "planned_runs": EXPECTED_RUNS,
            "verified_succeeded_runs": succeeded,
            "executed_succeeded_runs": counts["succeeded"],
            "skipped_verified_successes": counts["skipped_verified_success"],
            "failed_runs": counts["failed"],
            "not_started_runs": counts["not_started"],
            "workers": arguments.workers,
            "allocated_cpus": arguments.allocated_cpus,
            "allocated_memory_gb": arguments.allocated_memory_gb,
            "continuation_policy": CONTINUATION_POLICY,
            "sealed_features_opened": False,
            "sealed_labels_opened": False,
            "outcome_records_parsed_or_joined": False,
            "records": ordered,
        }
        write_json_exclusive(staging / "bundle_receipt.json", receipt)
        artifact_hash = freeze_tree(
            staging,
            {
                "artifact_class": "cobolt_v12_cpu_baseline_bundle_receipt",
                "campaign_id": CAMPAIGN_ID,
                "candidate_manifest_sha256": arguments.candidate_sha256,
                "stability_manifest_sha256": arguments.stability_sha256,
                "status": status,
                "planned_runs": EXPECTED_RUNS,
                "verified_succeeded_runs": succeeded,
                "outcomes_opened": False,
            },
        )
        publish_directory_noreplace(staging, target)
        print(
            json.dumps(
                {
                    "receipt": target.as_posix(),
                    "artifacts_sha256": artifact_hash,
                    "status": status,
                    "verified_succeeded_runs": succeeded,
                },
                sort_keys=True,
            )
        )
        return 0 if status == "succeeded" else 1
    finally:
        os.close(lock_descriptor)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--candidate-sha256", required=True)
    parser.add_argument("--stability", type=Path, required=True)
    parser.add_argument("--stability-sha256", required=True)
    parser.add_argument("--executor-sha256", required=True)
    parser.add_argument("--workers", type=int, default=MAX_WORKERS)
    parser.add_argument("--allocated-cpus", type=int, required=True)
    parser.add_argument("--allocated-memory-gb", type=int, required=True)
    parser.add_argument("--execution-root", type=Path)
    parser.add_argument("--receipt-root", type=Path)
    parser.add_argument("--execution-id")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--execute", action="store_true")
    return parser


def main() -> int:
    arguments = build_parser().parse_args()
    _assert_sha256(arguments.executor_sha256, "executor binding")
    if arguments.executor_sha256 != sha256_file(Path(__file__).resolve()):
        raise CpuBundleError("executor SHA-256 differs from the reviewed binding")
    if arguments.preflight:
        _, _, _, runs, startup = validate_startup(
            candidate_path=arguments.candidate,
            candidate_sha256=arguments.candidate_sha256,
            stability_path=arguments.stability,
            stability_sha256=arguments.stability_sha256,
            workers=arguments.workers,
            allocated_cpus=arguments.allocated_cpus,
            allocated_memory_gb=arguments.allocated_memory_gb,
        )
        print(
            json.dumps(
                {
                    **startup,
                    "runs": len(runs),
                    "continuation_policy": CONTINUATION_POLICY,
                    "execution_started": False,
                },
                sort_keys=True,
            )
        )
        return 0
    if (
        arguments.execution_root is None
        or arguments.receipt_root is None
        or arguments.execution_id is None
    ):
        raise CpuBundleError(
            "--execution-root, --receipt-root, and --execution-id are required with --execute"
        )
    return execute_bundle(arguments)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CpuBundleError as error:
        print(f"CPU bundle error: {error}", file=sys.stderr)
        raise SystemExit(2) from error
