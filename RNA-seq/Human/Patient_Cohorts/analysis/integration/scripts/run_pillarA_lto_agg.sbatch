#!/bin/bash
#SBATCH --job-name=pillarA-lto-agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarA-lto-agg_%j.out
#SBATCH --error=logs/pillarA-lto-agg_%j.err

# Pillar A LTO aggregation — combine leave-two-out results.

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Pillar A LTO aggregation on $(hostname)"
for LFC_VAL in 0 0.1 0.3 0.5; do
    if [ "${LFC_VAL}" = "0" ]; then TAG=nolfc;
    elif [ "${LFC_VAL}" = "0.1" ]; then TAG=lfc1;
    elif [ "${LFC_VAL}" = "0.3" ]; then TAG=lfc3;
    elif [ "${LFC_VAL}" = "0.5" ]; then TAG=lfc5;
    fi
    LFC_THR=${LFC_VAL} THR_TAG=${TAG} Rscript aggregate_pillar_A_lto.R
done
echo "[$(date)] Done LTO aggregation"
