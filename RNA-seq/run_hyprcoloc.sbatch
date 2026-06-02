#!/bin/bash
#SBATCH --job-name=hyprcoloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/hyprcoloc_%j.out
#SBATCH --error=RNA-seq/logs/hyprcoloc_%j.err

# ==============================================================================
# HyPrColoc multi-trait colocalization: Broadaway eQTL x ALL GWAS simultaneously
#
# Usage:
#   cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#   mkdir -p RNA-seq/logs
#   sbatch RNA-seq/run_hyprcoloc.sbatch
#
# Prerequisites:
#   - Phase 0D liftover: GWAS/MR_Data/hg19/*.tsv.gz must exist
#     (or at minimum PDFF native hg19 file)
#   - hyprcoloc installed in RNA-seq/.Rlib_causal/
#     If not: Rscript -e '.libPaths(c("RNA-seq/.Rlib_causal", .libPaths())); remotes::install_github("jrs95/hyprcoloc")'
# ==============================================================================

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

# Activate env BEFORE `set -u`: activate-binutils references an unbound ADDR2LINE under set -u
# (same workaround as run_208_subtype_coloc.sh). Previously failed in 0s with:
# "activate-binutils_linux-64.sh: line 68: ADDR2LINE: unbound variable".
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -euo pipefail

cd "${BASE}"

echo "=== HyPrColoc sbatch ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "Start: $(date)"
echo ""

Rscript RNA-seq/35u_hyprcoloc_broadaway.R

echo ""
echo "=== HyPrColoc complete ==="
echo "End: $(date)"
