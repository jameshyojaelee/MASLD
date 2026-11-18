#!/bin/bash
#SBATCH --job-name=c3_c3b_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/CONF1-20261109T221846Z/logs/c3_c3b_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/CONF1-20261109T221846Z/logs/c3_c3b_20261118_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CONF_DIR="${PROJECT_ROOT}/Analysis/Hypotheses/CONF1-20261109T221846Z"
mkdir -p "${CONF_DIR}/logs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  "${CONF_DIR}/scripts/c3_partition_20261118.R" \
  "${CONF_DIR}/out_c3_partition_20261118"

/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  "${CONF_DIR}/scripts/c3b_ambient_20261118.R" \
  "${CONF_DIR}/out_c3b_ambient_20261118"

echo "C3_AND_C3B_COMPLETE"
