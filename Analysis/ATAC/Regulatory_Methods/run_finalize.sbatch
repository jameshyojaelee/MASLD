#!/usr/bin/env bash
#SBATCH --job-name=masld_atac_finalize
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/finalize_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/finalize_%j.log
set -euo pipefail
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods
EXECUTION=${1:?Usage: sbatch run_finalize.sbatch EXECUTION_ROOT NEW_OUTPUT}
OUTPUT=${2:?New output required}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1
SOURCE=${REGULATORY_SOURCE_ROOT:-$TASK_ROOT/scripts}
/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python "$SOURCE/finalize_candidate.py" --execution-root "$EXECUTION" --output "$OUTPUT" --source-root "$SOURCE"
