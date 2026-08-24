"""Capture exact software and accelerator provenance for an immutable run."""

from __future__ import annotations

import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Sequence

from .artifacts import write_json_exclusive, write_text_exclusive


def _capture(command: Sequence[str]) -> tuple[int | None, str]:
    try:
        process = subprocess.run(
            list(command),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return None, f"command_not_found: {command[0]}\n"
    return process.returncode, process.stdout


def capture_runtime_lock(
    output_dir: str | Path,
    *,
    include_r: bool = False,
    include_cuda: bool = True,
) -> dict[str, object]:
    """Write provenance files before model code runs; unavailable tools are explicit."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    commands: dict[str, Sequence[str]] = {
        "pip-freeze.txt": (sys.executable, "-m", "pip", "freeze"),
    }
    if os.environ.get("CONDA_PREFIX"):
        commands["conda-explicit.txt"] = ("conda", "list", "--explicit")
    if include_cuda:
        commands["nvidia-smi.txt"] = ("nvidia-smi", "-q")
    if include_r:
        commands["R-sessionInfo.txt"] = (
            "Rscript",
            "--vanilla",
            "-e",
            "sessionInfo()",
        )
    command_status: dict[str, int | None] = {}
    for filename, command in commands.items():
        exit_code, output_text = _capture(command)
        command_status[filename] = exit_code
        write_text_exclusive(
            output / filename,
            output_text + f"\n# command={' '.join(command)}\n# exit_code={exit_code}\n",
        )
    metadata = {
        "schema_version": "masld-bench-runtime-lock-v1",
        "python_executable": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "hostname": platform.node(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_partition": os.environ.get("SLURM_JOB_PARTITION"),
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "command_status": command_status,
    }
    write_json_exclusive(output / "runtime.json", metadata)
    return metadata
