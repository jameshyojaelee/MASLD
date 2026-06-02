#!/bin/bash
#SBATCH --job-name=liver_xspecies_reanalysis
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/reanalysis_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/reanalysis_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00

# =============================================================================
# Targeted reanalysis: scripts 01 → 05 → 04
#
# Purpose: Regenerate gene concordance results with the new disease_vs_ctrl
# anchor (dvc_* columns), rebuild the unified atlas with cross_anchor_tier,
# and regenerate all plots including new plots 13 and 14.
#
# Skips: 02a/02b/02c (pathway), 03a/03b (WGCNA/TF) — unchanged, slow.
# Those pre-computed result files are read by 04 and 05 as normal.
# =============================================================================

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts"
LOG_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs"
mkdir -p "$LOG_DIR"

echo "=== Cross-Species Reanalysis (scripts 01 → 05 → 04) ==="
echo "Date  : $(date)"
echo "Node  : $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo ""

# -----------------------------------------------------------------------
# Step 1: Gene-level concordance with new disease_vs_ctrl anchor
# Adds dvc_category, dvc_n_concordant, dvc_diets_concordant columns to
# gene_concordance_per_gene.csv
# -----------------------------------------------------------------------
echo "================================================================"
echo "STEP 1: Gene-Level Concordance (01_corrected_gene_concordance.R)"
echo "================================================================"
Rscript "$SCRIPT_DIR/01_corrected_gene_concordance.R"
if [ $? -ne 0 ]; then echo "FAILED: step 1"; exit 1; fi
echo ""

# -----------------------------------------------------------------------
# Step 2: Unified atlas with cross_anchor_tier + translatability update
# Adds cross_anchor_tier: Dual_Conserved / NASH_Conserved / MASLD_Conserved
# Dual_Conserved genes receive +0.05 translatability score bonus
# -----------------------------------------------------------------------
echo "================================================================"
echo "STEP 2: Unified Concordance Atlas (05_unified_concordance_atlas.R)"
echo "================================================================"
Rscript "$SCRIPT_DIR/05_unified_concordance_atlas.R"
if [ $? -ne 0 ]; then echo "FAILED: step 2"; exit 1; fi
echo ""

# -----------------------------------------------------------------------
# Step 3: All visualizations including new plots 13 and 14
#   13_cross_anchor_tiers.pdf     — Dual/NASH/MASLD/Other gene counts
#   14_anchor_comparison_distribution.pdf — concordant diets per gene per anchor
# -----------------------------------------------------------------------
echo "================================================================"
echo "STEP 3: Visualizations (04_concordance_visualizations.R)"
echo "================================================================"
Rscript "$SCRIPT_DIR/04_concordance_visualizations.R"
if [ $? -ne 0 ]; then echo "FAILED: step 3"; exit 1; fi
echo ""

echo "=== Reanalysis complete: $(date) ==="
echo ""
echo "New/updated outputs:"
echo "  results/gene_concordance_per_gene.csv  (+ dvc_* columns)"
echo "  results/concordance_atlas_unified.csv  (+ cross_anchor_tier)"
echo ""
echo "Plots regenerated:"
ls /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/plots/*.pdf 2>/dev/null | xargs -I{} basename {}
