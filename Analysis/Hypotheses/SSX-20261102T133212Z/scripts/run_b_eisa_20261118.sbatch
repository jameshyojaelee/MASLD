#!/bin/bash
#SBATCH --job-name=ssx_eisa_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/SSX-20261102T133212Z/logs/ssx_eisa_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/SSX-20261102T133212Z/logs/ssx_eisa_20261118_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript \
  Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_b_eisa_20261118.R
ln -sfn ../../SSX-20261102T133212Z/out_eisa_20261118 \
  Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/bulk/out_eisa_20261118
