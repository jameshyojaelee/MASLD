#!/usr/bin/env bash
#SBATCH --job-name=g256-smoke
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579
SOURCE=${ROOT}/executions/gse256398-source-h5-21069407
CROSSWALK=${ROOT}/executions/gse256398-gencode49-crosswalk-21069439
QC=${ROOT}/executions/gse256398-qc-membership-21069669
SENSITIVITY=${ROOT}/executions/gse256398-qc-sensitivity-21069726
TARGET=${ROOT}/executions/gse256398-unlabeled-smoke-${SLURM_JOB_ID}

if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite GSE256398 smoke view: %s\n' "${TARGET}" >&2
  exit 1
fi
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH=${RUNTIME}/site-packages:${ROOT}/src:${ROOT}
python -m unittest -v tests.unit.test_build_gse256398_unlabeled_smoke
python "${ROOT}/scripts/build_gse256398_unlabeled_smoke.py" \
  --source "${SOURCE}" \
  --crosswalk "${CROSSWALK}" \
  --qc "${QC}" \
  --sensitivity "${SENSITIVITY}" \
  --output "${TARGET}"
python -m pip freeze > "${TARGET}/../gse256398-unlabeled-smoke-${SLURM_JOB_ID}.pip-freeze.txt"
sha256sum "${TARGET}/ARTIFACTS.json"
