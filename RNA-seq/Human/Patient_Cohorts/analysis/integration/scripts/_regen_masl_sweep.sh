#!/usr/bin/env bash
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/mash_loo_sweep

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
set -euo pipefail
cd "${SCRIPTS}"

export CONTRAST=masl_vs_healthy SLURM_CPUS_PER_TASK=4
Rscript lfc_sweep_loo.R \
  > "${LOGS}/regen_masl_sweep.out" 2> "${LOGS}/regen_masl_sweep.err"
