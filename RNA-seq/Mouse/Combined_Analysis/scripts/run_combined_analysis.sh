#!/bin/bash
# run_combined_analysis.sh
# Master script to run the combined MCD RNA-seq analysis pipeline

set -euo pipefail

# Use absolute path for SLURM compatibility
ROOT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Combined_Analysis"

# Use the unified mamba environment
MAMBA_EXE="/gpfs/commons/home/jameslee/.local/bin/micromamba"
ENV_PATH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/unified_env"

echo "=================================================="
echo "Combined MCD RNA-seq Analysis Pipeline"
echo "=================================================="
echo "Working directory: $ROOT_DIR"
echo "Environment: $ENV_PATH"
echo ""

cd "$ROOT_DIR"

echo "[Step 1/4] Preparing combined metadata..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/01_prepare_combined_metadata.R
echo ""

echo "[Step 2/4] Merging count matrices..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/02_merge_counts.R
echo ""

echo "[Step 3/4] Running DESeq2 with batch correction..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/03_run_deseq2_batch_corrected.R
echo ""

echo "[Step 4/4] Generating plots..."
$MAMBA_EXE run -p $ENV_PATH Rscript scripts/04_generate_plots.R
echo ""

echo "=================================================="
echo "Pipeline complete!"
echo "=================================================="
echo "Outputs:"
echo "  - metadata/samples.tsv"
echo "  - counts/combined_gene_counts.tsv"
echo "  - analysis_mcd_vs_control/deseq2_mcd_vs_control.tsv"
echo "  - analysis_mcd_vs_control/normalized_counts.csv"
echo "  - analysis_mcd_vs_control/vst_counts.csv"
echo "  - plots/*.pdf"
