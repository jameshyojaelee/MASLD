#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-diag-504
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
RUNTIME_PREFIX=/gpfs/commons/home/jameslee/micromamba/envs/immune-atlas
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SPLIT=${ROOT}/executions/model-data-052-21077595/split
TARGET=${ROOT}/executions/model-cpu-diag-504-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite elastic-net diagnostic artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}:${ROOT}/src" \
  python -m py_compile scripts/diagnose_elastic_net_convergence_study_50000.py
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}:${ROOT}/src" \
  OMP_NUM_THREADS="${OMP_NUM_THREADS}" MKL_NUM_THREADS="${MKL_NUM_THREADS}" \
  OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS}" \
  python scripts/diagnose_elastic_net_convergence_study_50000.py \
  --source "${SOURCE}" \
  --source-artifacts-sha256 672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509 \
  --split "${SPLIT}" \
  --split-artifacts-sha256 10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680 \
  --output "${STAGE}/diagnostic" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
sha256sum \
  scripts/diagnose_elastic_net_convergence_study_50000.py \
  slurm/diagnose_elastic_net_convergence_study_50000.sbatch \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
value = json.loads((root / "diagnostic/diagnostic.json").read_text())
if value.get("metrics_calculated") is not False or value.get("test_partition_predicted") is not False:
    raise SystemExit("diagnostic crossed the evaluation firewall")
freeze_tree(
    root,
    {
        "artifact_class": "cell_baseline_elastic_net_convergence_diagnostic",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "outer_fold": 0,
        "model_selected": False,
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
