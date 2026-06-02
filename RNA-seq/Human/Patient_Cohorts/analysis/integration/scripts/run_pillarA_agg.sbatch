#!/bin/bash
#SBATCH --job-name=pillarA-agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarA-agg_%j.out
#SBATCH --error=logs/pillarA-agg_%j.err

# Pillar A aggregation — combine CPSS + bootstrap + k-fold results.

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Pillar A aggregation on $(hostname)"

# Primary thresholds
for LFC_VAL in 0 0.1 0.3 0.5; do
    if [ "${LFC_VAL}" = "0" ]; then TAG=nolfc;
    elif [ "${LFC_VAL}" = "0.1" ]; then TAG=lfc1;
    elif [ "${LFC_VAL}" = "0.3" ]; then TAG=lfc3;
    elif [ "${LFC_VAL}" = "0.5" ]; then TAG=lfc5;
    fi
    echo "--- LFC=${LFC_VAL} (tag=${TAG}) ---"
    LFC_THR=${LFC_VAL} THR_TAG=${TAG} Rscript aggregate_pillar_A.R
done

echo "--- Extended LFC sweep ---"
LFC_THR=0.75 THR_TAG=lfc075 Rscript aggregate_pillar_A.R
LFC_THR=1.0  THR_TAG=lfc100 Rscript aggregate_pillar_A.R
LFC_THR=2.0  THR_TAG=lfc200 Rscript aggregate_pillar_A.R

echo "--- Threshold sweep comparison ---"
Rscript compare_lfc_threshold_sweep.R

echo "[$(date)] Done Pillar A aggregation"
