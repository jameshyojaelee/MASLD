#!/usr/bin/env bash
#SBATCH --job-name=model-data-914
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/pooled_graded_fibrosis_substrate-20260828T112808Z
mkdir -p "${OUT}"
${PYTHON} -m pip freeze > "${OUT}/pip_freeze.txt"
${PYTHON} -c "import sys,platform,json;print(json.dumps({'python':sys.version,'platform':platform.platform()},indent=2))" > "${OUT}/session_info.json"
${PYTHON} /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/scripts/substrate/freeze_pooled_graded_fibrosis_substrate.py "${OUT}"
echo "OUTPUT_DIR=${OUT}"
