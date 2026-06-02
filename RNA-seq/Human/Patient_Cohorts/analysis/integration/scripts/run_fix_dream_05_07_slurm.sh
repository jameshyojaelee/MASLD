#!/bin/bash
#SBATCH --job-name=liver_fix_dream
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fix_dream_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fix_dream_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=04:00:00

# =============================================================================
# Fix Dream Mega-Analysis and Restore Consensus DEGs
#
# Background: The Feb 25 6-cohort run (job 13908713) modified script 05 to
# include Hepatocytes + Macrophages as deconvolution covariates in the dream
# formula. These covariates are strongly correlated with disease status (MASLD
# causes hepatocyte loss + immune infiltration), so conditioning on them absorbed
# the disease signal — leaving only 66 DEGs and 0 Tier 1 consensus genes.
#
# Fix: Script 05 now uses the correct PRIMARY formula:
#       ~ group_binary + inferred_sex + (1|dataset)
# Deconvolution attribution analysis (adjusted vs unadjusted comparison) is
# handled separately by script 25_deconv_attribution.R.
#
# This script re-runs only:
#   05 (dream — now fixed) → 07 (consensus DEGs) → 08 (pathway) → 10 (volcano)
#   → 12 (library intersection)
# Scripts 00-04 and 06 do NOT need to be re-run (their outputs are correct).
# =============================================================================

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "============================================================"
echo "  DREAM FIX + CONSENSUS RESTORE (Scripts 05, 07, 08, 10, 12)"
echo "  Started: $(date)"
echo "  SLURM_JOB_ID: $SLURM_JOB_ID"
echo "  Formula fix: removed Hepatocytes + Macrophages from primary dream"
echo "============================================================"
echo ""

echo "=== 05: Dream Mega-Analysis (${SLURM_CPUS_PER_TASK} CPUs) ==="
echo "    Formula: ~ group_binary + inferred_sex + (1|dataset)"
Rscript analysis/integration/scripts/05_dream_mega_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05_dream_mega_analysis.R"; exit 1; fi

echo ""
echo "=== 07: Consensus DEGs ==="
Rscript analysis/integration/scripts/07_consensus_degs.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 07_consensus_degs.R"; exit 1; fi

echo ""
echo "=== 08: Pathway Analysis ==="
Rscript analysis/integration/scripts/08_pathway_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 08_pathway_analysis.R"; exit 1; fi

echo ""
echo "=== 10: Volcano Plots ==="
Rscript analysis/integration/scripts/10_volcano_plots.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 10_volcano_plots.R"; exit 1; fi

echo ""
echo "=== 12: Library Intersection ==="
Rscript analysis/integration/scripts/12_library_intersection.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 12_library_intersection.R"; exit 1; fi

echo ""
echo "============================================================"
echo "  DREAM FIX COMPLETE"
echo "  Finished: $(date)"
echo "  Expected: ~10K+ dream DEGs, Tier 1 > 100 genes"
echo "============================================================"
