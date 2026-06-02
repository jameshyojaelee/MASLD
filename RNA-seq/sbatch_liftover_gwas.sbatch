#!/bin/bash
#SBATCH --job-name=liftover_gwas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/liftover_gwas_%j.out
#SBATCH --error=RNA-seq/logs/liftover_gwas_%j.err

# =============================================================================
# Phase 0D: Liftover GWAS summary stats from hg38 to hg19
# =============================================================================

set -euo pipefail
echo "=== GWAS Liftover hg38 → hg19 ==="
echo "Start: $(date)"

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript RNA-seq/60_liftover_gwas_hg19.R

echo "=== Liftover complete ==="
echo "End: $(date)"
