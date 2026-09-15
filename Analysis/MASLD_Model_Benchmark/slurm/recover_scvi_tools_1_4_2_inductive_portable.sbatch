#!/usr/bin/env bash
#SBATCH --job-name=model-runtime-144
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
FAILED=${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-21079449.failed
FAILED_SOURCE_SHA256=b91862d9e82297ef9fe1b8cdc2036837280261d900a975f57c947191671d52ae
EXPECTED_SCVI_WHEEL_SHA256=c453b75b0fa0d222bf0eb3a531607d17f5ea0dbe5b848ab12fd8cc24599884e5
TARGET=${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-portable-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
RECOVERY_FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${RECOVERY_FAILED}" ]]; then
    mv "${STAGE}" "${RECOVERY_FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ ! -d "${FAILED}" || -e "${TARGET}" || -e "${STAGE}" || -e "${RECOVERY_FAILED}" ]]; then
  printf 'Portable scvi-tools 1.4.2 recovery path contract differs.\n' >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/downloads" \
  "${STAGE}/validation"

printf '%s  %s\n' "${FAILED_SOURCE_SHA256}" "${FAILED}/source.sha256" \
  | sha256sum --check --strict
cd "${ROOT}"
head -n 4 "${FAILED}/source.sha256" | sha256sum --check --strict
SCVI_WHEELS=("${FAILED}"/downloads/scvi_tools-1.4.2-*.whl)
if [[ "${#SCVI_WHEELS[@]}" -ne 1 || ! -f "${SCVI_WHEELS[0]}" ]]; then
  printf 'Expected exactly one recovered scvi-tools 1.4.2 wheel.\n' >&2
  exit 1
fi
printf '%s  %s\n' "${EXPECTED_SCVI_WHEEL_SHA256}" "${SCVI_WHEELS[0]}" \
  | sha256sum --check --strict
if find "${FAILED}/venv/lib/python3.12/site-packages" -type l -print -quit \
  | grep -q .; then
  printf 'Recovered scvi-tools 1.4.2 site-packages contains a symlink.\n' >&2
  exit 1
fi
cp -a "${FAILED}/venv/lib/python3.12/site-packages" "${STAGE}/site-packages"
cp -a "${SCVI_WHEELS[0]}" "${STAGE}/downloads/"
cp -a \
  "${FAILED}/validation/pip-check.txt" \
  "${FAILED}/validation/pip-freeze.txt" \
  "${FAILED}/validation/pip-install-report.json" \
  "${FAILED}/validation/pip-list.json" \
  "${FAILED}/validation/unit-tests.txt" \
  "${STAGE}/validation/"
cp -a \
  "${FAILED}/validation/inductive-probe" \
  "${STAGE}/validation/original-inductive-probe"

module purge
module load python/3.12.3-GCCcore-13.3.0
cd "${ROOT}"
export PYTHONPATH=${STAGE}/site-packages:${ROOT}/src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
python -m unittest tests.unit.test_scvi_scanvi_inductive \
  > "${STAGE}/validation/portable-unit-tests.txt" \
  2>&1
python scripts/probe_scvi_scanvi_1_4_2_inductive.py \
  --output "${STAGE}/validation/portable-inductive-probe" \
  > "${STAGE}/validation/portable-inductive-probe.stdout.json" \
  2> "${STAGE}/validation/portable-inductive-probe.stderr.txt"
python - "${STAGE}" <<'PY' > "${STAGE}/validation/portable-import-probe.json"
from pathlib import Path
import json
import platform
import sys

import anndata
import scvi
import sklearn
import torch

root = Path(sys.argv[1]).resolve(strict=True)
if scvi.__version__ != "1.4.2" or sys.version_info[:2] != (3, 12):
    raise SystemExit("portable scvi-tools or Python version differs")
print(json.dumps({
    "schema_version": "masld-bench-scvi-tools-portable-runtime-probe-v1",
    "status": "pass",
    "python": platform.python_version(),
    "scvi_tools": scvi.__version__,
    "torch": torch.__version__,
    "torch_cuda_build": torch.version.cuda,
    "cuda_visible_on_cpu_node": torch.cuda.is_available(),
    "anndata": anndata.__version__,
    "scikit_learn": sklearn.__version__,
    "portable_site_packages": True,
    "project_data_read": False,
    "sealed_outcomes_read": False,
}, sort_keys=True, separators=(",", ":")))
PY
sha256sum \
  slurm/recover_scvi_tools_1_4_2_inductive_portable.sbatch \
  scripts/probe_scvi_scanvi_1_4_2_inductive.py \
  src/masld_bench/adapters/scvi_scanvi_inductive.py \
  tests/unit/test_scvi_scanvi_inductive.py \
  "${FAILED}/source.sha256" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - \
  "${STAGE}" \
  "${EXPECTED_SCVI_WHEEL_SHA256}" \
  "${FAILED_SOURCE_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads(
    (root / "validation/portable-inductive-probe/probe.json").read_text(
        encoding="utf-8"
    )
)
if (
    probe.get("status") != "pass"
    or probe.get("scvi_tools") != "1.4.2"
    or probe.get("query_training_executed") is not False
    or probe.get("query_adaptation_executed") is not False
    or probe.get("scvi_state_unchanged") is not True
    or probe.get("scanvi_state_unchanged") is not True
    or probe.get("project_data_read") is not False
    or probe.get("sealed_outcomes_read") is not False
):
    raise SystemExit("portable scvi-tools 1.4.2 inductive probe differs")
freeze_tree(
    root / "validation/portable-inductive-probe",
    {
        "artifact_class": "scvi_scanvi_portable_inductive_query_probe",
        "scvi_tools_version": "1.4.2",
        "query_training_executed": False,
        "query_adaptation_executed": False,
        "learned_state_unchanged": True,
        "project_data_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
freeze_tree(
    root,
    {
        "artifact_class": "scvi_tools_1_4_2_inductive_portable_runtime",
        "model_roster": ["scvi_baseline", "scanvi_baseline"],
        "scvi_tools_version": "1.4.2",
        "scvi_tools_wheel_sha256": sys.argv[2],
        "code_commit": "dc95497",
        "code_license": "BSD-3-Clause",
        "recovered_failed_campaign_source_sha256": sys.argv[3],
        "failed_campaign_reason": "venv_launcher_symlinks_forbidden_by_immutable_tree",
        "portable_site_packages": True,
        "direct_inductive_query_probe": "passed",
        "query_training_allowed": False,
        "query_adaptation_allowed": False,
        "project_data_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Portable scvi-tools 1.4.2 inductive runtime: %s\n' "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/validation/portable-inductive-probe/ARTIFACTS.json"
