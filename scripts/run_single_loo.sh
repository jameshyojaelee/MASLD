#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --time=3:00:00

# Single LOO-CV dream job. HELD_OUT env var must be set.
set -euo pipefail

export PATH="/gpfs/commons/home/jameslee/miniforge3/condabin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/dream_loo_cv.R
