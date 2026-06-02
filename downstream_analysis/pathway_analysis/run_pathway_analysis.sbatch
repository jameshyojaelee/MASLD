#!/bin/bash
#SBATCH --job-name=pathway_analysis
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/Pathway/logs/pathway_analysis_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/Pathway/logs/pathway_analysis_%j.err

# ==============================================================================
# Pathway Analysis Pipeline - Master SLURM Script
# ==============================================================================
# This script orchestrates the complete pathway analysis workflow:
# 1. Environment setup (if needed)
# 2. Data preparation (extract log2FC values)
# 3. Gene set preparation (extract collections from MSigDB)
# 4. ORA analysis (Over-Representation Analysis)
# 5. GSEA analysis (Gene Set Enrichment Analysis)
# 6. Summary visualization generation
# ==============================================================================

set -e

# Configuration
PATHWAY_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/pathway_analysis"
SCRIPTS_DIR="${PATHWAY_DIR}/scripts"
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
ENV_PREFIX="${PATHWAY_DIR}/.mamba/pathway_analysis"

# Logging
LOG_DIR="${PATHWAY_DIR}/logs"
mkdir -p "${LOG_DIR}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/pipeline_${TIMESTAMP}.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

log "=========================================="
log "Pathway Analysis Pipeline Started"
log "=========================================="
log "Job ID: ${SLURM_JOB_ID:-local}"
log "Working directory: ${PATHWAY_DIR}"
log "Environment: ${ENV_PREFIX}"

# ==============================================================================
# Step 0: Check/Create Environment
# ==============================================================================
log ""
log "[0/5] Checking environment..."

if [ ! -d "${ENV_PREFIX}" ]; then
    log "Environment not found. Creating..."
    bash "${SCRIPTS_DIR}/environment_setup.sh" 2>&1 | tee -a "${LOG_FILE}"
else
    log "Environment exists at ${ENV_PREFIX}"
fi

# Set up environment for running R
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Create a wrapper for R to exclude user library
R_WRAPPER="${PATHWAY_DIR}/scripts/.r_wrapper.sh"
cat > "${R_WRAPPER}" << 'WRAPPER_EOF'
#!/bin/bash
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
ENV_PREFIX="/gpfs/commons/groups/sanjana_lab/Cas13/Pathway/.mamba/pathway_analysis"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
$MICROMAMBA run -p "${ENV_PREFIX}" Rscript --vanilla -e ".libPaths(.libPaths()[!grepl('jameslee', .libPaths())])" "$@"
WRAPPER_EOF
chmod +x "${R_WRAPPER}"

# ==============================================================================
# Step 1: Prepare Data
# ==============================================================================
log ""
log "[1/5] Preparing data (extracting log2FC values)..."

${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/prepare_data.R" 2>&1 | tee -a "${LOG_FILE}"

if [ $? -eq 0 ]; then
    log "Data preparation completed successfully"
else
    log "ERROR: Data preparation failed"
    exit 1
fi

# ==============================================================================
# Step 2: Prepare Gene Sets
# ==============================================================================
log ""
log "[2/5] Preparing gene sets from MSigDB..."

${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/prepare_genesets.R" 2>&1 | tee -a "${LOG_FILE}"

if [ $? -eq 0 ]; then
    log "Gene set preparation completed successfully"
else
    log "ERROR: Gene set preparation failed"
    exit 1
fi

# ==============================================================================
# Step 3: Run ORA Analysis
# ==============================================================================
log ""
log "[3/5] Running Over-Representation Analysis (ORA)..."

${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/run_ora.R" 2>&1 | tee -a "${LOG_FILE}"

if [ $? -eq 0 ]; then
    log "ORA analysis completed successfully"
else
    log "ERROR: ORA analysis failed"
    exit 1
fi

# ==============================================================================
# Step 4: Run GSEA Analysis
# ==============================================================================
log ""
log "[4/5] Running Gene Set Enrichment Analysis (GSEA)..."

${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/run_gsea.R" 2>&1 | tee -a "${LOG_FILE}"

if [ $? -eq 0 ]; then
    log "GSEA analysis completed successfully"
else
    log "ERROR: GSEA analysis failed"
    exit 1
fi

# ==============================================================================
# Step 5: Generate Summary Plots
# ==============================================================================
log ""
log "[5/5] Generating summary visualizations..."

${MICROMAMBA} run -p "${ENV_PREFIX}" Rscript "${SCRIPTS_DIR}/generate_summary_plots.R" 2>&1 | tee -a "${LOG_FILE}"

if [ $? -eq 0 ]; then
    log "Summary visualization completed successfully"
else
    log "ERROR: Summary visualization failed"
    exit 1
fi

# ==============================================================================
# Summary
# ==============================================================================
log ""
log "=========================================="
log "Pipeline Completed Successfully!"
log "=========================================="
log ""
log "Results locations:"
log "  ORA results:     ${PATHWAY_DIR}/results/ora/"
log "  GSEA results:    ${PATHWAY_DIR}/results/gsea/"
log "  Individual plots: ${PATHWAY_DIR}/plots/individual/"
log "  Summary plots:    ${PATHWAY_DIR}/plots/summary/"
log ""
log "Key output files:"
log "  - all_ora_results.RData"
log "  - all_gsea_results.RData"
log "  - pathway_summary_multipanel_ora.pdf"
log "  - gsea_nes_heatmap.pdf"
log "  - pathway_counts_comparison.pdf"
log ""
log "Pipeline log: ${LOG_FILE}"
