#!/bin/bash
# run_meta_analysis.sh
# Run all meta-analysis scripts sequentially

set -e

MAMBA_EXE="/gpfs/commons/home/jameslee/.local/bin/micromamba"
ENV_PATH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/unified_env"
WORK_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Combined_Analysis"

cd "$WORK_DIR"

echo "=================================================="
echo "Combined MCD Meta-Analysis Pipeline"
echo "=================================================="
echo "Working directory: $WORK_DIR"
echo "Environment: $ENV_PATH"
echo ""

# Step 1: Count DEGs
echo "[Step 1/5] Counting upregulated DEGs across datasets..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/meta_01_count_degs.R
echo ""

# Step 2: Overlap analysis
echo "[Step 2/5] Computing overlap metrics..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/meta_02_overlap_analysis.R
echo ""

# Step 3: LFC correlation
echo "[Step 3/5] Generating LFC correlation plots..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/meta_03_lfc_correlation.R
echo ""

# Step 4: Human comparison
echo "[Step 4/5] Cross-species comparison..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/meta_04_human_comparison.R
echo ""

# Step 5: Library comparison
echo "[Step 5/5] Library coverage analysis..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/meta_05_library_comparison.R
echo ""

echo "=================================================="
echo "Meta-Analysis Complete!"
echo "=================================================="
echo "Outputs in: $WORK_DIR/meta_analysis/"
echo "Plots in: $WORK_DIR/plots/"
