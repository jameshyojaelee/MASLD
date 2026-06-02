#!/bin/bash
#SBATCH --job-name=pillarC-agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarC-agg_%j.out
#SBATCH --error=logs/pillarC-agg_%j.err

# Pillar C aggregation — combine permutation null results.

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Pillar C aggregation on $(hostname)"
for LFC_VAL in 0 0.1 0.3 0.5; do
    if [ "${LFC_VAL}" = "0" ]; then TAG=nolfc;
    elif [ "${LFC_VAL}" = "0.1" ]; then TAG=lfc1;
    elif [ "${LFC_VAL}" = "0.3" ]; then TAG=lfc3;
    elif [ "${LFC_VAL}" = "0.5" ]; then TAG=lfc5;
    fi
    echo "--- LFC=${LFC_VAL} (tag=${TAG}) ---"
    LFC_THR=${LFC_VAL} THR_TAG=${TAG} Rscript aggregate_pillar_C.R
done
echo "[$(date)] Done Pillar C aggregation"
