#!/usr/bin/env bash
#SBATCH --job-name=pseudobulk_de
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=Analysis/SingleCell/integration/logs/pseudobulk_de_%j.out
#SBATCH --error=Analysis/SingleCell/integration/logs/pseudobulk_de_%j.err

set -euo pipefail

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Pseudobulk DE (limma-voom) ==="
echo "Date: $(date)"

Rscript Analysis/SingleCell/scripts/pseudobulk_de.R

echo ""
echo "=== Pseudobulk DE Complete ==="
echo "End: $(date)"
