#!/bin/bash
# 05: GSE281367 pseudobulk (19, snapatac2) -> pooled n=30 DA (18, rnaseq), atomic.
#SBATCH --job-name=edgeR
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/pbda_%j.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/pbda_%j.err
set -eo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

echo "[$(date)] === STEP 19: GSE281367 pseudobulk (snapatac2) ==="
micromamba activate snapatac2
python Analysis/ATAC/Human_External/scripts/19_pseudobulk_gse281367.py

echo "[$(date)] === STEP 18: pooled n=30 DA (rnaseq) ==="
micromamba deactivate
micromamba activate rnaseq
Rscript Analysis/ATAC/Human_External/scripts/18_pooled_da_gse281367.R
echo "[$(date)] === DONE ==="
