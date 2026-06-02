#!/bin/bash
#SBATCH --job-name=figS_varpart
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figS_varpart_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figS_varpart_%j.err

set -o pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
mkdir -p scripts/figures/logs

source "$HOME/.bashrc"
micromamba activate rnaseq
set -u

echo "=== figS_variance_partition: regen at canonical FIGS_SENS_DIR ==="
echo "Started: $(date)"

Rscript scripts/figures/figS_variance_partition.R

echo "=== DONE at $(date) ==="
