#!/usr/bin/env python3
"""Freeze one existing micromamba environment without installing packages."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Sequence

from masld_bench.artifacts import (
    canonical_hash,
    freeze_tree,
    publish_directory_noreplace,
    sha256_file,
    write_json_exclusive,
    write_text_exclusive,
)


class EnvironmentFreezeError(RuntimeError):
    """Raised when an environment cannot be frozen exactly."""


_PROBE = r'''
import importlib.metadata
import json
import os
import platform
import sys

def normalized(name):
    return name.strip().lower().replace("_", "-")

records = sorted({
    f"{normalized(distribution.metadata.get('Name', ''))}=={distribution.version}"
    for distribution in importlib.metadata.distributions()
    if distribution.metadata.get('Name')
})
critical_names = (
    "anndata",
    "h5py",
    "numpy",
    "pandas",
    "scanpy",
    "scikit-learn",
    "scipy",
)
critical = {
    name: importlib.metadata.version(name)
    for name in critical_names
}
print(json.dumps({
    "critical_versions": critical,
    "distribution_records": records,
    "python_executable": os.path.realpath(sys.executable),
    "python_implementation": platform.python_implementation(),
    "python_prefix": os.path.realpath(sys.prefix),
    "python_version": platform.python_version(),
    "user_site_enabled": bool(__import__('site').ENABLE_USER_SITE),
}, sort_keys=True, separators=(",", ":")))
'''


def _run(command: Sequence[str], *, label: str) -> str:
    environment = dict(os.environ)
    environment["PYTHONNOUSERSITE"] = "1"
    process = subprocess.run(
        list(command),
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=environment,
    )
    if process.returncode != 0:
        raise EnvironmentFreezeError(
            f"{label} failed with exit {process.returncode}: {process.stderr[-4000:]}"
        )
    return process.stdout


def _conda_meta_inventory(prefix: Path) -> list[dict[str, Any]]:
    metadata = prefix / "conda-meta"
    if not metadata.is_dir():
        raise EnvironmentFreezeError(f"environment has no conda-meta directory: {prefix}")
    records: list[dict[str, Any]] = []
    for path in sorted(metadata.glob("*.json")):
        if path.is_symlink() or not path.is_file():
            raise EnvironmentFreezeError(f"invalid conda metadata member: {path}")
        records.append(
            {
                "path": path.name,
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )
    if not records:
        raise EnvironmentFreezeError(f"environment has no conda package records: {prefix}")
    return records


def freeze_environment(
    *,
    micromamba: Path,
    prefix: Path,
    runtime_id: str,
    output_root: Path,
) -> Path:
    micromamba = micromamba.resolve(strict=True)
    prefix = prefix.resolve(strict=True)
    if not micromamba.is_file() or not os.access(micromamba, os.X_OK):
        raise EnvironmentFreezeError(f"micromamba is not executable: {micromamba}")
    if not prefix.is_dir():
        raise EnvironmentFreezeError(f"environment prefix is not a directory: {prefix}")
    output_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{runtime_id}.", dir=output_root))
    try:
        explicit = _run(
            (
                micromamba.as_posix(),
                "list",
                "--explicit",
                "--sha256",
                "--no-rc",
                "--prefix",
                prefix.as_posix(),
            ),
            label="micromamba explicit export",
        )
        pip_freeze = _run(
            (
                micromamba.as_posix(),
                "run",
                "--prefix",
                prefix.as_posix(),
                "env",
                "PYTHONNOUSERSITE=1",
                "python",
                "-m",
                "pip",
                "freeze",
                "--all",
            ),
            label="pip freeze",
        )
        probe_text = _run(
            (
                micromamba.as_posix(),
                "run",
                "--prefix",
                prefix.as_posix(),
                "env",
                "PYTHONNOUSERSITE=1",
                "python",
                "-c",
                _PROBE,
            ),
            label="isolated Python probe",
        ).strip()
        try:
            probe = json.loads(probe_text)
        except json.JSONDecodeError as error:
            raise EnvironmentFreezeError("isolated Python probe returned invalid JSON") from error
        if not isinstance(probe, dict) or probe.get("user_site_enabled") is not False:
            raise EnvironmentFreezeError("isolated Python probe did not disable user site")
        if probe.get("python_prefix") != prefix.as_posix():
            raise EnvironmentFreezeError("isolated Python prefix differs from requested prefix")
        distribution_records = probe.pop("distribution_records", None)
        if not isinstance(distribution_records, list) or not distribution_records:
            raise EnvironmentFreezeError("isolated Python distribution inventory is empty")
        conda_meta = _conda_meta_inventory(prefix)

        explicit_path = staging / "conda-explicit.txt"
        pip_path = staging / "pip-freeze.txt"
        distribution_path = staging / "python-distributions.txt"
        conda_meta_path = staging / "conda-meta-inventory.json"
        write_text_exclusive(explicit_path, explicit.rstrip() + "\n")
        write_text_exclusive(pip_path, pip_freeze.rstrip() + "\n")
        write_text_exclusive(
            distribution_path,
            "\n".join(str(item) for item in distribution_records) + "\n",
        )
        write_json_exclusive(conda_meta_path, conda_meta)

        lock = {
            "schema_version": "masld-bench-python-environment-lock-v1",
            "runtime_id": runtime_id,
            "environment_prefix": prefix.as_posix(),
            "micromamba_path": micromamba.as_posix(),
            "micromamba_sha256": sha256_file(micromamba),
            "micromamba_version": _run(
                (micromamba.as_posix(), "--version"), label="micromamba version"
            ).strip(),
            "python_executable": probe["python_executable"],
            "python_implementation": probe["python_implementation"],
            "python_version": probe["python_version"],
            "python_no_user_site_required": True,
            "critical_versions": probe["critical_versions"],
            "conda_explicit_sha256": sha256_file(explicit_path),
            "pip_freeze_sha256": sha256_file(pip_path),
            "python_distributions_sha256": sha256_file(distribution_path),
            "conda_meta_inventory_sha256": sha256_file(conda_meta_path),
            "conda_meta_content_sha256": canonical_hash(conda_meta),
            "package_installation_performed": False,
        }
        lock_path = staging / "environment.lock.json"
        write_json_exclusive(lock_path, lock)
        lock_sha256 = sha256_file(lock_path)
        target = output_root / f"{runtime_id}--{lock_sha256[:16]}"
        freeze_tree(
            staging,
            {
                "artifact_class": "python_environment_lock",
                "runtime_id": runtime_id,
                "environment_lock_sha256": lock_sha256,
            },
        )
        publish_directory_noreplace(staging, target)
        return target
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--micromamba", required=True, type=Path)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--runtime-id", required=True)
    parser.add_argument("--output-root", required=True, type=Path)
    arguments = parser.parse_args(argv)
    target = freeze_environment(
        micromamba=arguments.micromamba,
        prefix=arguments.prefix,
        runtime_id=arguments.runtime_id,
        output_root=arguments.output_root,
    )
    lock = target / "environment.lock.json"
    print(
        json.dumps(
            {
                "environment_artifact": {
                    "path": lock.as_posix(),
                    "sha256": sha256_file(lock),
                    "size_bytes": lock.stat().st_size,
                    "media_type": "application/json",
                    "role": f"environment:{arguments.runtime_id}",
                },
                "lock_directory": target.as_posix(),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
