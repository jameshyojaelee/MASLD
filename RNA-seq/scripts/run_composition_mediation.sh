#!/bin/bash
#SBATCH --partition=cpu --qos=interactive
#SBATCH --cpus-per-task=16 --mem=128G --time=48:00:00
#SBATCH --job-name=med_F_trans
#SBATCH --output=RNA-seq/logs/med_F_trans_%j.out
#SBATCH --error=RNA-seq/logs/med_F_trans_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export MED_N_TOP=${MED_N_TOP:-200}
export MED_N_BOOT=${MED_N_BOOT:-1000}

Rscript RNA-seq/scripts/run_composition_mediation.R
