#!/usr/bin/env bash
#SBATCH --job-name=masld-perturb-runtime
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --array=0-6
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT="${PROJECT}/Analysis/MASLD_Model_Benchmark"
LANES=(gears scgpt scgen cpa cellot genepert perturblib)
ENVS=(perturbation_gears perturbation_scgpt perturbation_scgen_cpa perturbation_scgen_cpa perturbation_gears perturbation_state perturbation_state)
LANE=${LANES[${SLURM_ARRAY_TASK_ID}]}
ENV_NAME=${ENVS[${SLURM_ARRAY_TASK_ID}]}
PYTHON="${PROJECT}/.envs/${ENV_NAME}/bin/python"
TARGET="${ROOT}/executions/perturbation-runtime-${LANE}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
STAGE="${TARGET}.staging"
if [[ ! -x "${PYTHON}" || -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Runtime audit precondition differs: %s\n' "${LANE}" >&2
  exit 1
fi
cd "${ROOT}"
"${PYTHON}" scripts/audit_existing_perturbation_runtime.py \
  --lane "${LANE}" \
  --output "${STAGE}"
sha256sum \
  scripts/audit_existing_perturbation_runtime.py \
  slurm/audit_existing_perturbation_runtimes.sbatch \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${STAGE}" "${LANE}" "${ENV_NAME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "existing_perturbation_runtime_audit",
    "model_lane": sys.argv[2],
    "environment": sys.argv[3],
    "installation_performed": False,
    "forward_fixture_performed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
