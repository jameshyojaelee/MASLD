#!/usr/bin/env bash
#SBATCH --job-name=model-probe-015
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --gres=gpu:b6k:1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Can the scVI/scANVI code path run under ENFORCED determinism on sm_120?
#
# Job 21108408 observed bit-identical cross-process reproducibility, but the probe
# does not enable torch.use_deterministic_algorithms, so that reproducibility is a
# property of the run, not a requirement. Production needs results reproducible BY
# REQUIREMENTS. This probe answers whether the requirements are available.
#
# It applies the same knobs Cobolt uses -- use_deterministic_algorithms(True),
# cudnn disabled/deterministic/benchmark-off, CUBLAS_WORKSPACE_CONFIG=:4096:8 --
# and then runs the SHA-bound probe's own run() entry point unmodified.
#
# If some op has no deterministic kernel on sm_120, torch raises naming that op.
# This probe CAPTURES the op name rather than falling back silently, so the
# decision is explicit. A failure here is a finding about scvi-tools 1.4.2 on
# sm_120, NOT about B6K.
#
# For contrast it also records the warn_only=True outcome, which reports the same
# ops as warnings instead of raising, so a single run enumerates every offending
# op rather than stopping at the first.
#
# Synthetic 40-gene / 120-row fixture only. No project data, no development
# outcome, no held-back outcome, no benchmark metric. Direct sbatch; never enters the
# dispatcher roster. Inherited job 21100293 untouched.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-portable-21079606"
RUNTIME_ARTIFACTS_SHA256=d885ada24fdad9f6e3d57a92ff870ca8ff05f07701ac6b0762f55a8cb19db61a
PROBE="${ROOT}/scripts/probe_scvi_scanvi_1_4_2_inductive.py"
PROBE_SHA256=dccd8ae73483ae3c7e0047a01fd1d1301b42ccba9a11ebac149814b3d7df4362

TARGET="${ROOT}/executions/scvi-scanvi-determinism-feasibility-b6k-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite determinism probe artifact: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}"

cd "${ROOT}"
sha256sum --check --strict <<EOF
${RUNTIME_ARTIFACTS_SHA256}  ${RUNTIME}/ARTIFACTS.json
${PROBE_SHA256}  ${PROBE}
EOF

nvidia-smi > "${STAGE}/nvidia-smi.txt" 2>&1 || true

module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH=${RUNTIME}/site-packages:${ROOT}:${ROOT}/src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

python - "${STAGE}" "${PROBE}" <<'PY'
from __future__ import annotations

import importlib.util
import json
import os
import re
import sys
import traceback
import warnings
from pathlib import Path

stage = Path(sys.argv[1]).resolve(strict=True)
probe_path = Path(sys.argv[2]).resolve(strict=True)

import torch

result: dict[str, object] = {
    "schema_version": "masld-bench-scvi-scanvi-determinism-feasibility-v1",
    "question": (
        "Can the scVI/scANVI code path run under ENFORCED determinism on sm_120, so "
        "production results are reproducible by contract rather than by observation?"
    ),
    "attribution_rule": (
        "A failure here is a finding about the scvi-tools 1.4.2 code path on sm_120, "
        "NOT about B6K hardware."
    ),
    "torch_version": torch.__version__,
    "cuda_runtime": torch.version.cuda,
    "cublas_workspace_config": os.environ.get("CUBLAS_WORKSPACE_CONFIG"),
}
if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
    raise SystemExit("exactly one visible CUDA device is required")
result["device_name"] = torch.cuda.get_device_name(0)
result["device_capability"] = list(torch.cuda.get_device_capability(0))
result["arch_list"] = list(torch.cuda.get_arch_list())

spec = importlib.util.spec_from_file_location("scvi_scanvi_probe_under_determinism", probe_path)
probe = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = probe
spec.loader.exec_module(probe)

_OP = re.compile(r"([A-Za-z0-9_.]+)\s+does not have a deterministic implementation")


def apply_knobs(warn_only: bool) -> None:
    torch.use_deterministic_algorithms(True, warn_only=warn_only)
    torch.backends.cudnn.enabled = False
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


# --- Pass 1: strict. Does anything raise? ----------------------------------
strict: dict[str, object] = {"warn_only": False}
try:
    apply_knobs(warn_only=False)
    probe.run(stage / "strict-probe")
    strict["status"] = "pass"
    strict["nondeterministic_ops"] = []
except Exception as error:  # noqa: BLE001
    text = f"{type(error).__name__}: {error}"
    strict["status"] = "fail"
    strict["error"] = text[:2000]
    strict["nondeterministic_ops"] = sorted(set(_OP.findall(text)))
    (stage / "strict-traceback.txt").write_text(
        "".join(traceback.format_exception(type(error), error, error.__traceback__)),
        encoding="utf-8",
    )
result["strict"] = strict

# --- Pass 2: warn_only. Enumerate EVERY offending op, not just the first ----
warn: dict[str, object] = {"warn_only": True}
try:
    apply_knobs(warn_only=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        probe.run(stage / "warn-only-probe")
    messages = [str(w.message) for w in caught]
    ops = sorted({op for m in messages for op in _OP.findall(m)})
    warn["status"] = "pass"
    warn["nondeterministic_ops"] = ops
    warn["warning_count"] = len(messages)
    (stage / "warn-only-warnings.txt").write_text("\n".join(messages) + "\n", encoding="utf-8")
except Exception as error:  # noqa: BLE001
    warn["status"] = "fail"
    warn["error"] = f"{type(error).__name__}: {error}"[:2000]
result["warn_only"] = warn

# --- Recommendation, derived not asserted ----------------------------------
all_ops = sorted(set(strict.get("nondeterministic_ops") or []) | set(warn.get("nondeterministic_ops") or []))
result["nondeterministic_ops_union"] = all_ops
if strict.get("status") == "pass" and not all_ops:
    result["enforceable"] = True
    result["recommendation"] = (
        "ENABLE enforced determinism for production. The scVI/scANVI path runs clean "
        "under use_deterministic_algorithms(True) with cudnn disabled and "
        "CUBLAS_WORKSPACE_CONFIG=:4096:8 on sm_120. Reproducibility becomes a contract."
    )
else:
    result["enforceable"] = False
    result["recommendation"] = (
        "DO NOT silently fall back. The ops listed in nondeterministic_ops_union have "
        "no deterministic kernel on sm_120, so enforced determinism is unavailable for "
        "this code path. The choice -- accept observed-but-unenforced determinism for "
        "this lane, recorded as such -- is a decision to be taken explicitly, not a default."
    )
result["status"] = "pass"
(stage / "determinism.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2, sort_keys=True))
PY

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" "${PROBE_SHA256}" <<'PYC'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "scvi_scanvi_determinism_feasibility_b6k",
        "determinism_sha256": sha256_file(root / "determinism.json"),
        "probe_source_sha256": sys.argv[3],
        "runtime_artifacts_sha256": sys.argv[2],
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PYC

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'scvi/scanvi determinism feasibility complete: %s\n' "${TARGET}"
