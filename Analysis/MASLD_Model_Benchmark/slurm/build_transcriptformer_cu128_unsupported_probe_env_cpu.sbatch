#!/usr/bin/env bash
#SBATCH --job-name=model-probe-020
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# UNSUPPORTED CONFIGURATION. This environment exists only to measure how much
# TranscriptFormer's own outputs move between torch 2.5.1 and torch 2.8.0 with
# the GPU pool held fixed. It is not a production environment and nothing else
# may build on it.
#
# TranscriptFormer 0.6.1 pins torch==2.5.1 by exact equality, along with
# pytorch-lightning, numpy, anndata, scanpy, scipy, pandas and h5py. There is
# therefore no supported way to run it on torch 2.8.0+cu128. This build installs
# the package with --no-deps and substitutes torch deliberately, so the pin is
# violated on purpose and the violation is recorded rather than hidden.
#
# pip check is run WITHOUT masking its exit status and its full output is kept.
# A clean pip check here would mean the substitution did not happen.
#
# Reads no project data, no development outcome, no held-back outcome, and
# calculates no benchmark metric.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
micromamba="/gpfs/commons/home/jameslee/.local/bin/micromamba"
reference_env="${root}/executions/environments/transcriptformer-torch2.5.1-21070738"
wheel_source="${reference_env}/wheels/transcriptformer-0.6.1-py3-none-any.whl"
wheel_sha256=b3fb64ebc238b5316fa4c564c63186172acf6c3d2a914a4bfbf46452311c8cdf

target="${root}/executions/environments/transcriptformer-cu128-UNSUPPORTED-PROBE-${SLURM_JOB_ID}"
stage="${target}.staging"
environment="${stage}/env"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}/validation" "${stage}/wheels"

# Reuse the already-built wheel rather than rebuilding from source, so the
# package itself is bit-identical to the one the supported environment uses and
# the only difference between the two environments is the torch substitution.
printf '%s  %s\n' "${wheel_sha256}" "${wheel_source}" | sha256sum --check --strict
cp "${wheel_source}" "${stage}/wheels/"
wheel="${stage}/wheels/transcriptformer-0.6.1-py3-none-any.whl"

"${micromamba}" create --yes --prefix "${environment}" \
  --channel conda-forge \
  python=3.11 \
  pip

python_bin="${environment}/bin/python"

# torch first, from the cu128 index, so nothing later can pull a different build.
"${python_bin}" -m pip install --no-cache-dir \
  --index-url https://download.pytorch.org/whl/cu128 \
  torch==2.8.0

# The rest of TranscriptFormer's declared pins, at their declared versions,
# with torch deliberately absent from this list.
"${python_bin}" -m pip install --no-cache-dir \
  anndata==0.11.4 \
  h5py==3.14 \
  hydra-core==1.3.2 \
  numpy==2.2.6 \
  pandas==2.2.2 \
  psutil==5.9 \
  scanpy==1.11.2 \
  scipy==1.15.3

# pytorch-lightning declares a torch range rather than an exact pin, so it is
# installed without dependencies to guarantee it cannot pull torch 2.5.1 back in.
"${python_bin}" -m pip install --no-cache-dir --no-deps pytorch-lightning==2.5.1
"${python_bin}" -m pip install --no-cache-dir \
  lightning-utilities torchmetrics

# The remaining declared pins, so the only requirement pip check reports is the
# torch substitution itself rather than a list of absences that mean nothing.
# cellxgene-census is deliberately NOT installed: it exists for remote corpus
# access, local inference from a local checkpoint never calls it, and pulling it
# would drag a large dependency set that could perturb the pinned stack we are
# trying to hold fixed. Its absence is expected in the pip check output.
"${python_bin}" -m pip install --no-cache-dir \
  boto3==1.38.27 \
  pynvml==12 \
  pytest==8.4 \
  requests==2.32.3 \
  timeout-decorator==0.5

# The package itself, deps suppressed. This is the deliberate violation.
"${python_bin}" -m pip install --no-cache-dir --no-deps "${wheel}"

# Record the broken requirements verbatim. pip check exits non-zero when it
# finds one, which is the expected outcome here, so its status is captured
# rather than allowed to end the job.
pip_check_status=0
"${python_bin}" -m pip check > "${stage}/validation/pip-check.txt" 2>&1 \
  || pip_check_status=$?
printf 'pip check exit status: %s\n' "${pip_check_status}" \
  >> "${stage}/validation/pip-check.txt"
"${python_bin}" -m pip freeze > "${stage}/validation/pip-freeze.txt"

"${python_bin}" - "${stage}" <<'PY' > "${stage}/validation/import-probe.json"
import json
import sys
from pathlib import Path

import torch

record = {
    "torch_version": torch.__version__,
    "torch_cuda_build": torch.version.cuda,
    "torch_arch_list": list(torch.cuda.get_arch_list()),
    "cuda_available_on_cpu_build_node": torch.cuda.is_available(),
    "configuration": "UNSUPPORTED",
    "declared_torch_pin": "2.5.1",
    "installed_torch_differs_from_pin": torch.__version__.split("+")[0] != "2.5.1",
    "production_use_permitted": False,
    "project_data_read": False,
    "sealed_outcomes_read": False,
}
try:
    import transcriptformer

    record["transcriptformer_import"] = "pass"
    record["transcriptformer_version"] = getattr(
        transcriptformer, "__version__", "unknown"
    )
except Exception as error:
    record["transcriptformer_import"] = "fail"
    record["transcriptformer_import_error_type"] = type(error).__name__
    record["transcriptformer_import_error"] = str(error)[:800]
print(json.dumps(record, indent=2, sort_keys=True))
PY

"${micromamba}" list --explicit --sha256 --no-rc --prefix "${environment}" \
  > "${stage}/conda-explicit.txt"
sha256sum "${stage}"/wheels/*.whl > "${stage}/validation/wheels.sha256"

cat > "${stage}/PROBE_ONLY.json" <<EOF
{
  "schema_version": "masld-bench-unsupported-configuration-probe-env-v1",
  "status": "unsupported_configuration_probe_environment",
  "production_use_permitted": false,
  "may_produce_benchmark_metrics": false,
  "purpose": "Measure how far TranscriptFormer's own outputs move between torch 2.5.1 and torch 2.8.0 with the GPU pool held fixed. Nothing else may build on this environment.",
  "declared_pin_violated": "transcriptformer 0.6.1 requires torch==2.5.1",
  "installed_torch": "2.8.0+cu128",
  "wheel_sha256": "${wheel_sha256}",
  "reference_supported_environment": "executions/environments/transcriptformer-torch2.5.1-21070738",
  "pip_check_exit_status": ${pip_check_status},
  "authorised_by": "team lead, 2026-08-30, probe only, explicitly not the production re-extraction"
}
EOF

# Digest the recorded files individually rather than freezing the whole tree.
# A conda prefix contains symlinks (env/bin/2to3 among others) and the immutable
# tree helper refuses symlinks outright, which is why the supported build writes
# a SHA256SUMS over named files instead. Same approach here, same reason.
(
  cd "${stage}"
  sha256sum \
    PROBE_ONLY.json \
    conda-explicit.txt \
    validation/import-probe.json \
    validation/pip-check.txt \
    validation/pip-freeze.txt \
    validation/wheels.sha256 \
    > SHA256SUMS
)

mv "${stage}" "${target}"
printf 'built unsupported-configuration probe environment: %s\n' "${target}"
cat "${target}/SHA256SUMS"
