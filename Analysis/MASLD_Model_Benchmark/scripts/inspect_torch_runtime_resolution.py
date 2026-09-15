#!/usr/bin/env python3
"""Report which torch and which CUDA libraries a stacked PYTHONPATH resolves to.

The campaign runs its pinned ``torch 2.6.0+cu124`` overlay on top of a base
layer that also ships a torch, and the overlay's ``effective-pip-freeze.txt``
lists both CUDA 12 and CUDA 13 NVIDIA packages.  Static ELF metadata says the
CUDA 13 entries are a freeze-time union of the stack rather than something the
process can load, because torch's libraries carry an ``$ORIGIN``-relative RPATH
into their own layer and a ``libcudart.so.12`` soname.  This script confirms
that inside a real interpreter instead of arguing it from packaging rules.

It reads no biological data, no outcome, and computes no benchmark metric.
"""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
from typing import Any


SCHEMA = "masld-bench-torch-runtime-resolution-v1"
CUDA_LIBRARY = re.compile(
    r"lib(?:cudart|cudnn[a-z_]*|cublas[A-Za-z]*|cufft|curand|cusolver|cusparse"
    r"|nccl|nvrtc|nvJitLink|torch[a-z_]*|c10[a-z_]*)\.so[.0-9]*$"
)
# The CUDA 13 distributions whose presence in the overlay freeze prompted this
# check.  They are looked up by name, not assumed absent.
CUDA13_DISTRIBUTIONS = (
    "nvidia-cuda-runtime",
    "nvidia-cudnn-cu13",
    "nvidia-cublas",
    "nvidia-nccl-cu13",
    "nvidia-cusparselt-cu13",
)
CUDA12_DISTRIBUTIONS = (
    "nvidia-cuda-runtime-cu12",
    "nvidia-cudnn-cu12",
    "nvidia-cublas-cu12",
    "nvidia-nccl-cu12",
)


def readelf_dynamic(path: Path) -> dict[str, Any]:
    """DT_NEEDED / RPATH / RUNPATH for one shared object."""

    try:
        result = subprocess.run(
            ["readelf", "-d", str(path)],
            check=True,
            text=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        return {"status": f"unavailable: {type(error).__name__}"}
    needed: list[str] = []
    rpath: list[str] = []
    runpath: list[str] = []
    for line in result.stdout.splitlines():
        match = re.search(r"\(([A-Z]+)\).*\[(.*)\]", line)
        if match is None:
            continue
        tag, value = match.group(1), match.group(2)
        if tag == "NEEDED":
            needed.append(value)
        elif tag == "RPATH":
            rpath.extend(value.split(":"))
        elif tag == "RUNPATH":
            runpath.extend(value.split(":"))
    return {
        "status": "read",
        "needed": needed,
        "rpath": rpath,
        "runpath": runpath,
        # RPATH takes precedence over LD_LIBRARY_PATH; RUNPATH does not.  Which
        # one is present decides whether the environment can redirect these.
        "search_is_overridable_by_ld_library_path": not rpath and bool(runpath),
    }


def loaded_cuda_libraries() -> list[str]:
    """Shared objects actually mapped into this process, from /proc/self/maps."""

    try:
        text = Path("/proc/self/maps").read_text(encoding="utf-8")
    except OSError:
        return []
    paths = {
        line.rsplit(" ", 1)[-1].strip()
        for line in text.splitlines()
        if " /" in line
    }
    return sorted(path for path in paths if CUDA_LIBRARY.search(path))


def distribution_report(names: tuple[str, ...]) -> dict[str, Any]:
    report: dict[str, Any] = {}
    for name in names:
        try:
            distribution = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            report[name] = {"found": False}
            continue
        report[name] = {
            "found": True,
            "version": distribution.version,
            "location": str(getattr(distribution, "_path", "unknown")),
        }
    return report


def torch_copies_on_path() -> list[dict[str, str]]:
    """Every importable torch on sys.path, in resolution order.

    A second torch further down the path is shadowed only by element order.  If
    that order ever changes, jobs would silently switch runtime, so the copies
    are enumerated rather than assumed unique.
    """

    copies: list[dict[str, str]] = []
    for entry in sys.path:
        candidate = Path(entry) / "torch" / "version.py"
        if not candidate.is_file():
            continue
        text = candidate.read_text(encoding="utf-8")
        version = re.search(r"__version__ = '([^']+)'", text)
        cuda = re.search(r"cuda: Optional\[str\] = '([^']+)'", text)
        copies.append(
            {
                "sys_path_entry": entry,
                "version": version.group(1) if version else "unknown",
                "cuda": cuda.group(1) if cuda else "none",
            }
        )
    return copies


def checkpoint_report(path: Path) -> dict[str, Any]:
    import torch

    record: dict[str, Any] = {"path": str(path)}
    try:
        state = torch.load(path, map_location="cpu", weights_only=True)
    except Exception as error:  # noqa: BLE001 - the failure is the measurement
        record["weights_only_load"] = "failed"
        record["error_type"] = type(error).__name__
        record["error"] = str(error)
        return record
    record["weights_only_load"] = "passed"
    if not isinstance(state, dict):
        record["payload_type"] = type(state).__name__
        return record
    tensors = {k: v for k, v in state.items() if isinstance(v, torch.Tensor)}
    record["entry_count"] = len(state)
    record["tensor_count"] = len(tensors)
    record["all_entries_are_tensors"] = len(tensors) == len(state)
    record["non_tensor_keys"] = sorted(set(state) - set(tensors))
    record["parameter_count"] = int(sum(v.numel() for v in tensors.values()))
    record["dtypes"] = sorted({str(v.dtype) for v in tensors.values()})
    record["keys"] = sorted(tensors)[:64]
    record["key_count"] = len(tensors)
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--checkpoint", action="append", type=Path, default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    # sys.path is captured before importing torch so the record shows the
    # resolution order the interpreter was actually given.
    resolution_order = list(sys.path)

    import torch

    torch_lib = Path(torch.__file__).resolve().parent / "lib"
    report: dict[str, Any] = {
        "schema_version": SCHEMA,
        "label": args.label,
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "node": platform.node(),
        "sys_path": resolution_order,
        "torch_copies_on_sys_path": torch_copies_on_path(),
        "resolved_torch": {
            "version": torch.__version__,
            "file": str(Path(torch.__file__).resolve()),
            "cuda_build": torch.version.cuda,
        },
        "shadowed_torch_count": max(0, len(torch_copies_on_path()) - 1),
        "nvidia_namespace_path": [],
        "cuda13_distributions": distribution_report(CUDA13_DISTRIBUTIONS),
        "cuda12_distributions": distribution_report(CUDA12_DISTRIBUTIONS),
        "libtorch_global_deps": readelf_dynamic(torch_lib / "libtorch_global_deps.so"),
        "libtorch_cuda": readelf_dynamic(torch_lib / "libtorch_cuda.so"),
        "loaded_cuda_libraries_after_import": loaded_cuda_libraries(),
        "outcomes_read": False,
        "project_data_read": False,
        "metrics_calculated": False,
    }
    try:
        report["torch_arch_flags"] = torch._C._cuda_getArchFlags().split()
    except Exception as error:  # noqa: BLE001
        report["torch_arch_flags"] = f"unavailable: {type(error).__name__}"
    try:
        import nvidia

        report["nvidia_namespace_path"] = list(getattr(nvidia, "__path__", []))
    except ImportError:
        report["nvidia_namespace_path"] = []

    report["checkpoints"] = [checkpoint_report(path) for path in args.checkpoint]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
