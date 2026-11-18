#!/bin/bash
#SBATCH --job-name=eqtl_psd_C
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/compare_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/compare_%j.err
set -eo pipefail
eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export RUN_DIR=${MASLD_PROJECT_ROOT}/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118
Rscript ${RUN_DIR}/code/compare.R
