#!/bin/bash
#SBATCH --job-name=snapatac2
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/coembed_%j.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/coembed_%j.err
set -eo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"; export MASLD_PROJECT_ROOT="$PROJ"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
echo "[$(date)] co-embed GSE244832 + GSE281367 on $(hostname)"
python Analysis/ATAC/Human_External/scripts/30_coembed_umap.py
echo "[$(date)] done"
