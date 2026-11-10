#!/bin/bash
#SBATCH --job-name=u8a_null
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive/logs/%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
U=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive
SEED=$1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla $U/02_null.R $SEED 19999 $U/out/null_seed_${SEED} 16
