#!/bin/bash
#SBATCH --job-name=pillarA_extlfc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pillarA_extlfc_%j.out
#SBATCH --error=logs/pillarA_extlfc_%j.err

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs

export MAMBA_EXE='/gpfs/commons/home/jameslee/.local/bin/micromamba'
export MAMBA_ROOT_PREFIX='/gpfs/commons/home/jameslee/micromamba'
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate rnaseq
set -eo pipefail

echo "=== LFC=0.75 (tag=lfc075) ==="
LFC_THR=0.75 THR_TAG=lfc075 Rscript aggregate_pillar_A.R

echo "=== LFC=1.0 (tag=lfc100) ==="
LFC_THR=1.0  THR_TAG=lfc100 Rscript aggregate_pillar_A.R

echo "=== LFC=2.0 (tag=lfc200) ==="
LFC_THR=2.0  THR_TAG=lfc200 Rscript aggregate_pillar_A.R

echo "Done."
