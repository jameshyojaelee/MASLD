#!/bin/bash
#SBATCH --job-name=rev_U0U1
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U0_U1/rederive/logs/%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U0_U1/rederive
OUT0=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U0_U1/rederive/out/u0
OUT1=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U0_U1/rederive/out/u1
/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/python -I /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U0_U1/rederive/u0_roster_recount.py $OUT0
/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/pip freeze > $OUT0/pip_freeze.txt 2>/dev/null || true
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U0_U1/rederive/u1_refit.R $OUT1 $OUT0
