#!/bin/bash
#SBATCH --job-name=atac_diffbind
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Mouse_Bulk/logs/diffbind_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Mouse_Bulk/logs/diffbind_%j.err

# =============================================================================
# Module 1 Step 2: DiffBind Differential Accessibility + ChIPseeker Annotation
# Runs after Snakemake pipeline completes (alignment → peaks → consensus)
# =============================================================================

set -euo pipefail

PIPELINE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Mouse_Bulk"
cd "$PIPELINE_DIR"

echo "============================================"
echo "DiffBind Differential Accessibility Analysis"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "============================================"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export PIPELINE_DIR="$PIPELINE_DIR"
export MC_CORES=4
Rscript scripts/03_diffbind_analysis.R

echo ""
echo "============================================"
echo "DiffBind analysis completed: $(date)"
echo "============================================"
