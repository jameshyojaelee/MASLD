#!/bin/bash
#SBATCH --job-name=s3a_loco_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3a_loco_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3a_loco_20261118_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "${PROJECT_ROOT}/Analysis/Hypotheses/SS-20261025T135237Z/logs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  "${PROJECT_ROOT}/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_loco_20261118.R" \
  --out "${PROJECT_ROOT}/Analysis/Hypotheses/SS-20261025T135237Z/out_s3a_loco_20261118"
