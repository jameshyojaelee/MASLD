#!/usr/bin/env bash
#SBATCH --job-name=tf-runtime
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PIP_DISABLE_PIP_VERSION_CHECK=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
source_root="${root}/executions/transcriptformer-source-21070659"
source_artifacts_sha256="d0ee2c7d5d9bef3441449e425bd920e206de9cd846af438c16e84d27e04b2d96"
micromamba="/gpfs/commons/home/jameslee/.local/bin/micromamba"
target="${root}/executions/environments/transcriptformer-torch2.5.1-${SLURM_JOB_ID}"
environment="${target}/env"
export MAMBA_ROOT_PREFIX="${target}/mamba-root"

if [[ -e "${target}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer runtime: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${target}/validation" "${target}/wheels"
printf '%s  %s\n' "${source_artifacts_sha256}" "${source_root}/ARTIFACTS.json" \
  | sha256sum --check --strict

"${micromamba}" create --yes --prefix "${environment}" \
  --channel conda-forge \
  python=3.11 \
  pip

"${environment}/bin/python" -m pip wheel \
  --no-deps \
  --wheel-dir "${target}/wheels" \
  "${source_root}/source"
wheel="$(find "${target}/wheels" -maxdepth 1 -type f -name 'transcriptformer-0.6.1-*.whl' -print -quit)"
if [[ -z "${wheel}" ]]; then
  printf 'Pinned source did not build the expected TranscriptFormer 0.6.1 wheel.\n' >&2
  exit 1
fi
"${environment}/bin/python" -m pip install --no-cache-dir "${wheel}"
"${environment}/bin/python" -m pip check \
  > "${target}/validation/pip-check.txt"

"${environment}/bin/python" - <<'PY' > "${target}/validation/import-probe.json"
from importlib.metadata import version
import json
import sys

import anndata
import h5py
import hydra
import numpy
import pandas
import pytorch_lightning
import scanpy
import scipy
import torch
import transcriptformer

expected = {
    "anndata": "0.11.4",
    "boto3": "1.38.27",
    "cellxgene-census": "1.17.0",
    "h5py": "3.14.0",
    "hydra-core": "1.3.2",
    "numpy": "2.2.6",
    "pandas": "2.2.2",
    "psutil": "5.9.0",
    "pynvml": "12.0.0",
    "pytest": "8.4.0",
    "pytorch-lightning": "2.5.1",
    "requests": "2.32.3",
    "scanpy": "1.11.2",
    "scipy": "1.15.3",
    "timeout-decorator": "0.5.0",
    "torch": "2.5.1",
    "transcriptformer": "0.6.1",
}
observed = {name: version(name) for name in expected}
if observed != expected:
    raise SystemExit(f"TranscriptFormer runtime versions differ: {observed!r}")
print(json.dumps({
    "cuda_available_on_cpu_preflight": bool(torch.cuda.is_available()),
    "python": sys.version,
    "status": "pass_cpu_import",
    "torch_cuda_build": torch.version.cuda,
    "versions": observed,
}, sort_keys=True, separators=(",", ":")))
PY

"${micromamba}" list --explicit --sha256 --no-rc --prefix "${environment}" \
  > "${target}/conda-explicit.txt"
"${environment}/bin/python" -m pip freeze --all \
  > "${target}/pip-freeze.txt"
sha256sum "${target}"/wheels/*.whl > "${target}/validation/wheels.sha256"
sha256sum \
  "${source_root}/ARTIFACTS.json" \
  "${root}/slurm/build_transcriptformer_runtime.sbatch" \
  > "${target}/source-manifest.sha256"

wheel_sha256="$(sha256sum "${wheel}" | cut -d' ' -f1)"
conda_sha256="$(sha256sum "${target}/conda-explicit.txt" | cut -d' ' -f1)"
pip_sha256="$(sha256sum "${target}/pip-freeze.txt" | cut -d' ' -f1)"
probe_sha256="$(sha256sum "${target}/validation/import-probe.json" | cut -d' ' -f1)"
cat > "${target}/ARTIFACTS.json" <<EOF
{
  "code_commit": "c943a89a4de9511a8ab1010715bf1cfa8827cc99",
  "conda_explicit_sha256": "${conda_sha256}",
  "environment_prefix": "${environment}",
  "import_probe_sha256": "${probe_sha256}",
  "pip_freeze_sha256": "${pip_sha256}",
  "python_no_user_site_required": true,
  "schema_version": "masld-bench-transcriptformer-runtime-v1",
  "status": "pass_cpu_runtime",
  "transcriptformer_version": "0.6.1",
  "wheel_sha256": "${wheel_sha256}"
}
EOF
(
  cd "${target}"
  sha256sum \
    ARTIFACTS.json \
    conda-explicit.txt \
    pip-freeze.txt \
    source-manifest.sha256 \
    validation/import-probe.json \
    validation/pip-check.txt \
    validation/wheels.sha256 \
    > SHA256SUMS
)
chmod -R a-w "${target}"
printf 'TranscriptFormer isolated runtime candidate: %s\n' "${target}"
