#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=348_lr_bulk_concord
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/348_lr_bulk_concord_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/348_lr_bulk_concord_%j.err

set -eo pipefail

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "[start] $(date)"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/348_stage_ccc_bulk_concordance.R
echo "[done] $(date)"
