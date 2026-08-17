"""Create and execute hash-locked GPU run specifications."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path, sha256_tree


ALLOWED_GPU_SCRIPTS = {
    "train_reference.py", "train_reference_sensitivity.py", "train_de_novo.py",
    "train_update.py", "train_sequence.py",
}
PATH_ARGUMENTS = {
    "--config", "--prepared", "--prepared-lock", "--contract-lock",
    "--reference-model", "--reference-embedding", "--selection-lock", "--refit-lock", "--output",
    "--bi-replay-unlock", "--policy",
}
RESULT_MANIFESTS = {
    "train_reference.py": "reference_manifest.json",
    "train_reference_sensitivity.py": "reference_manifest.json",
    "train_de_novo.py": "denovo_manifest.json",
    "train_update.py": "update_manifest.json",
    "train_sequence.py": "sequence_manifest.json",
}


def _canonical_digest(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def source_identity(pipeline_root: str | Path) -> dict[str, Any]:
    root = Path(pipeline_root).resolve()
    files: list[Path] = []
    for relative in ("masld_cl", "scripts", "tests", "slurm"):
        base = root / relative
        files.extend(
            path for path in base.rglob("*")
            if path.is_file()
            and "__pycache__" not in path.parts
            and path.suffix in {".py", ".sbatch"}
        )
    files.extend([root / "config_v1.json", root / "README.md"])
    files.extend(path for path in (root / "reference").rglob("*") if path.is_file())
    entries = []
    for path in sorted(set(files), key=lambda value: value.relative_to(root).as_posix()):
        entries.append({
            "path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_path(path),
        })
    if not entries:
        raise ContractError("GPU source identity contains no files")
    return {"entries": entries, "identity_sha256": _canonical_digest(entries)}


def verify_source_identity_payload(identity: dict[str, Any]) -> dict[str, Any]:
    """Verify a historical source identity without requiring it to equal current source."""
    if (
        not isinstance(identity, dict) or not identity.get("entries")
        or identity.get("identity_sha256") != _canonical_digest(identity.get("entries"))
    ):
        raise ContractError("historical source identity payload is invalid")
    return identity


def _normalize_arguments(arguments: list[str], project_root: Path) -> tuple[list[str], dict[str, str]]:
    normalized = list(arguments)
    paths: dict[str, str] = {}
    index = 0
    while index < len(normalized):
        flag = normalized[index]
        if flag in PATH_ARGUMENTS:
            if flag in paths or index + 1 >= len(normalized):
                raise ContractError(f"execution argument is missing or duplicated: {flag}")
            path = Path(normalized[index + 1])
            if not path.is_absolute():
                path = project_root / path
            path = path.resolve()
            normalized[index + 1] = str(path)
            paths[flag] = str(path)
            index += 2
        else:
            index += 1
    for required in ("--config", "--output"):
        if required not in paths:
            raise ContractError(f"GPU execution spec lacks {required}")
    return normalized, paths


def _artifact_identity(path: str | Path) -> dict[str, Any]:
    value = Path(path).resolve()
    if value.is_file():
        return {
            "kind": "file", "realpath": str(value),
            "size_bytes": value.stat().st_size, "sha256": sha256_path(value),
        }
    if value.is_dir():
        digest, files = sha256_tree(value)
        return {
            "kind": "directory", "realpath": str(value),
            "tree_sha256": digest, "files": files,
        }
    raise ContractError(f"locked execution input does not exist: {value}")


def execution_lock_payload(
    spec: dict[str, Any], spec_path: str | Path, pipeline_root: str | Path,
) -> dict[str, Any]:
    pipeline_root = Path(pipeline_root).resolve()
    project_root = pipeline_root.parents[2]
    arguments, paths = _normalize_arguments(spec["arguments"], project_root)
    inputs = {
        flag: _artifact_identity(path)
        for flag, path in sorted(paths.items()) if flag != "--output"
    }
    config_path = Path(paths["--config"])
    with config_path.open() as handle:
        config_payload = json.load(handle)
    return {
        "schema_version": "masld-cl-execution-lock-v2",
        "spec_realpath": str(Path(spec_path).resolve()),
        "spec_file_sha256": sha256_path(spec_path),
        "spec_sha256": spec["spec_sha256"],
        "normalized_arguments": arguments,
        "output_realpath": paths["--output"],
        "config_file_sha256": sha256_path(config_path),
        "config_content_sha256": _canonical_digest(config_payload),
        "source_identity": source_identity(pipeline_root),
        "input_artifacts": inputs,
    }


def write_execution_lock(
    spec_path: str | Path, pipeline_root: str | Path, output: str | Path,
) -> dict[str, Any]:
    with Path(spec_path).open() as handle:
        spec = json.load(handle)
    _validate_spec_content(spec)
    lock = execution_lock_payload(spec, spec_path, pipeline_root)
    lock["lock_sha256"] = _canonical_digest(lock)
    write_json_exclusive(output, lock)
    return lock


def _validate_spec_content(spec: dict[str, Any]) -> None:
    if spec.get("schema_version") not in {
        "masld-cl-execution-spec-v1", "masld-cl-execution-spec-v2"
    }:
        raise ContractError("unsupported execution spec")
    payload = {key: value for key, value in spec.items() if key != "spec_sha256"}
    if _canonical_digest(payload) != spec.get("spec_sha256"):
        raise ContractError("execution spec content hash mismatch")
    if spec.get("script") not in ALLOWED_GPU_SCRIPTS:
        raise ContractError(f"script is not allowed in a GPU execution spec: {spec.get('script')}")
    arguments = spec.get("arguments")
    if not isinstance(arguments, list) or not all(isinstance(x, str) for x in arguments):
        raise ContractError("execution arguments must be a JSON string array")


def _load_and_verify_execution_lock(
    spec: dict[str, Any], spec_path: Path, pipeline_root: Path,
) -> dict[str, Any]:
    lock_path = spec_path.with_suffix(spec_path.suffix + ".source-lock.json")
    if not lock_path.is_file():
        raise ContractError(f"GPU execution spec lacks pre-run source lock: {lock_path}")
    with lock_path.open() as handle:
        lock = json.load(handle)
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if _canonical_digest(payload) != lock.get("lock_sha256"):
        raise ContractError("execution source lock content hash mismatch")
    observed = execution_lock_payload(spec, spec_path, pipeline_root)
    for key in (
        "spec_realpath", "spec_file_sha256", "spec_sha256", "normalized_arguments",
        "output_realpath", "config_file_sha256", "config_content_sha256",
        "source_identity", "input_artifacts",
    ):
        if lock.get(key) != observed.get(key):
            raise ContractError(f"execution source or input changed after locking: {key}")
    return lock


def _load_and_verify_historical_execution_lock(
    spec: dict[str, Any], spec_path: Path, pipeline_root: Path,
) -> dict[str, Any]:
    """Verify a completed run's frozen lock without requiring old source as current."""
    lock_path = spec_path.with_suffix(spec_path.suffix + ".source-lock.json")
    if not lock_path.is_file():
        raise ContractError(f"GPU execution spec lacks pre-run source lock: {lock_path}")
    with lock_path.open() as handle:
        lock = json.load(handle)
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if _canonical_digest(payload) != lock.get("lock_sha256"):
        raise ContractError("execution source lock content hash mismatch")

    project_root = Path(pipeline_root).resolve().parents[2]
    arguments, paths = _normalize_arguments(spec["arguments"], project_root)
    config_path = Path(paths["--config"])
    with config_path.open() as handle:
        config_payload = json.load(handle)
    observed_static = {
        "spec_realpath": str(spec_path.resolve()),
        "spec_file_sha256": sha256_path(spec_path),
        "spec_sha256": spec["spec_sha256"],
        "normalized_arguments": arguments,
        "output_realpath": paths["--output"],
        "config_file_sha256": sha256_path(config_path),
        "config_content_sha256": _canonical_digest(config_payload),
        "input_artifacts": {
            flag: _artifact_identity(path)
            for flag, path in sorted(paths.items()) if flag != "--output"
        },
    }
    for key, value in observed_static.items():
        if lock.get(key) != value:
            raise ContractError(f"historical execution input changed after locking: {key}")
    identity = lock.get("source_identity") or {}
    if not identity.get("entries") or identity.get("identity_sha256") != _canonical_digest(
        identity.get("entries", [])
    ):
        raise ContractError("historical execution source identity is invalid")
    return lock


def execute_run_spec(
    spec_path: str | Path, pipeline_root: str | Path, record_dir: str | Path,
) -> dict[str, Any]:
    spec_path = Path(spec_path).resolve()
    pipeline_root = Path(pipeline_root).resolve()
    with spec_path.open() as handle:
        spec = json.load(handle)
    _validate_spec_content(spec)
    resource = os.environ.get("MASLD_RESOURCE_CLASS")
    if spec.get("resource_class") != resource:
        raise ContractError(
            f"execution spec requires {spec.get('resource_class')}, allocated {resource}"
        )
    lock = _load_and_verify_execution_lock(spec, spec_path, pipeline_root)
    script_name = spec["script"]
    script = pipeline_root / "scripts" / script_name
    output = Path(lock["output_realpath"])
    if output.exists():
        raise ContractError(f"GPU output path already exists: {output}")
    python = Path(os.environ["MASLD_PYTHON_BIN"]).resolve()
    record_dir = Path(record_dir)
    record_dir.mkdir(parents=True, exist_ok=True)
    if (record_dir / "execution_record.json").exists():
        raise ContractError("execution record directory was already used")
    started = datetime.now(timezone.utc).isoformat()
    command = [str(python), str(script), *lock["normalized_arguments"]]
    status = "failed"
    exit_code = None
    result_identity = None
    try:
        completed = subprocess.run(command, check=False)
        exit_code = completed.returncode
        if exit_code != 0:
            raise subprocess.CalledProcessError(exit_code, command)
        result_manifest = output / RESULT_MANIFESTS[script_name]
        if not result_manifest.is_file():
            raise ContractError(f"GPU run omitted result manifest: {result_manifest}")
        with result_manifest.open() as handle:
            result = json.load(handle)
        if result.get("config_sha256") != lock["config_content_sha256"]:
            raise ContractError("GPU result manifest config hash mismatch")
        output_hash, output_files = sha256_tree(output)
        result_identity = {
            "output_realpath": str(output.resolve()),
            "output_tree_sha256": output_hash,
            "output_files": output_files,
            "result_manifest_realpath": str(result_manifest.resolve()),
            "result_manifest_sha256": sha256_path(result_manifest),
        }
        status = "complete"
    finally:
        record = {
            "schema_version": "masld-cl-execution-record-v2",
            "config_sha256": lock["config_content_sha256"],
            "spec_realpath": str(spec_path),
            "spec_file_sha256": sha256_path(spec_path),
            "spec_sha256": spec["spec_sha256"],
            "execution_lock_realpath": str(
                spec_path.with_suffix(spec_path.suffix + ".source-lock.json").resolve()
            ),
            "execution_lock_sha256": sha256_path(
                spec_path.with_suffix(spec_path.suffix + ".source-lock.json")
            ),
            "source_identity_sha256": lock["source_identity"]["identity_sha256"],
            "input_artifacts": lock["input_artifacts"],
            "resource_class": resource,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_cpus": os.environ.get("SLURM_CPUS_PER_TASK"),
            "slurm_gpus": os.environ.get("SLURM_GPUS"),
            "started_utc": started,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "hostname": platform.node(),
            "status": status,
            "exit_code": exit_code,
            "command": command,
            "result_identity": result_identity,
        }
        write_json_exclusive(record_dir / "execution_record.json", record)
        (record_dir / ("EXECUTION_COMPLETE" if status == "complete" else "FAILED")).touch(
            exist_ok=False
        )
    return record


def verify_execution_record(
    record_path: str | Path, pipeline_root: str | Path, config_sha256: str,
    *, allow_historical_source: bool = False,
) -> dict[str, Any]:
    path = Path(record_path).resolve()
    with path.open() as handle:
        record = json.load(handle)
    if (
        record.get("schema_version") != "masld-cl-execution-record-v2"
        or record.get("status") != "complete"
        or record.get("exit_code") != 0
        or record.get("config_sha256") != config_sha256
        or not path.parent.joinpath("SLURM_COMPLETE").is_file()
        or path.parent.joinpath("FAILED.exit_code").exists()
    ):
        raise ContractError(f"GPU execution record is incomplete or failed: {path}")
    spec_path = Path(record["spec_realpath"])
    with spec_path.open() as handle:
        spec = json.load(handle)
    _validate_spec_content(spec)
    lock_path = Path(record["execution_lock_realpath"])
    if sha256_path(lock_path) != record.get("execution_lock_sha256"):
        raise ContractError("execution lock changed after GPU run")
    loader = (
        _load_and_verify_historical_execution_lock
        if allow_historical_source else _load_and_verify_execution_lock
    )
    lock = loader(spec, spec_path.resolve(), Path(pipeline_root).resolve())
    if lock["source_identity"]["identity_sha256"] != record.get("source_identity_sha256"):
        raise ContractError("GPU execution source identity differs from its record")
    result = record.get("result_identity") or {}
    output = Path(result.get("output_realpath", ""))
    tree_hash, files = sha256_tree(output)
    if tree_hash != result.get("output_tree_sha256") or files != result.get("output_files"):
        raise ContractError("GPU output bundle changed after execution")
    manifest = Path(result.get("result_manifest_realpath", ""))
    if sha256_path(manifest) != result.get("result_manifest_sha256"):
        raise ContractError("GPU result manifest changed after execution")
    return record


def execution_owned_files(record: dict[str, Any]) -> set[str]:
    """Return every absolute output path protected by a verified record."""
    identity = record.get("result_identity") or {}
    root = Path(identity.get("output_realpath", "")).resolve()
    files = {
        str((root / item["path"]).resolve())
        for item in identity.get("output_files", [])
    }
    manifest = identity.get("result_manifest_realpath")
    if manifest:
        files.add(str(Path(manifest).resolve()))
    return files


def require_execution_ownership(
    record: dict[str, Any], paths: list[str | Path], *, role: str,
) -> None:
    owned = execution_owned_files(record)
    missing = sorted(
        str(Path(path).resolve()) for path in paths
        if str(Path(path).resolve()) not in owned
    )
    if missing:
        raise ContractError(f"{role} artifacts are not owned by its GPU record: {missing}")
