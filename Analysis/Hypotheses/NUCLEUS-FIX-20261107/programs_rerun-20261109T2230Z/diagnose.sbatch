#!/bin/bash
#SBATCH --job-name=unique-prog-diag
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/NUCLEUS-FIX-20261107/programs_rerun-20261109T2230Z/diagnose_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript Analysis/Hypotheses/NUCLEUS-FIX-20261107/programs_rerun-20261109T2230Z/diagnose_two.R
