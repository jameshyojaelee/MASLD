#!/bin/bash
set -e

# Configuration
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DOWNSTREAM_DIR="${PROJECT_ROOT}/downstream_analysis"
PATHWAY_ENV="${DOWNSTREAM_DIR}/pathway_analysis/.mamba/pathway_analysis"
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"

# Argument Parsing
MODE="all"
if [[ "$#" -gt 0 ]]; then
    MODE="$1"
fi

# Help
if [[ "$MODE" == "--help" || "$MODE" == "-h" ]]; then
    echo "Usage: ./run_all_downstream.sh [OPTION]"
    echo ""
    echo "Options:"
    echo "  all           Run all analyses (default)"
    echo "  essentiality  Run only Essentiality Analysis"
    echo "  pathway       Run only Pathway Analysis"
    echo "  patient       Run only Patient Comparisons"
    echo ""
    exit 0
fi

echo "=========================================================="
echo "   MASLD Downstream Analysis Pipeline"
echo "   Mode: ${MODE}"
echo "=========================================================="
echo "Date: $(date)"
echo "Root: ${PROJECT_ROOT}"

# Check for Core DEGs
if [ ! -f "${PROJECT_ROOT}/final_core_degs.csv" ]; then
    echo "ERROR: final_core_degs.csv not found in project root."
    exit 1
fi
echo "[OK] Found final_core_degs.csv"

# ------------------------------------------------------------------
# 1. Essentiality Analysis
# ------------------------------------------------------------------
if [[ "$MODE" == "all" || "$MODE" == "essentiality" ]]; then
    echo ""
    echo "[1/3] Running Essentiality Analysis..."
    cd "${PROJECT_ROOT}"

    # Run master pipeline
    python3 Analysis/downstream_analysis/essentiality/scripts/run_essentiality_pipeline.py
    # Run supplementary plots
    python3 Analysis/downstream_analysis/essentiality/scripts/analyze_essentiality.py
    python3 Analysis/downstream_analysis/essentiality/scripts/extended_essentiality_analysis.py

    echo "[OK] Essentiality Analysis Complete."
else
    echo "[SKIP] Essentiality Analysis"
fi

# ------------------------------------------------------------------
# 2. Pathway Analysis
# ------------------------------------------------------------------
if [[ "$MODE" == "all" || "$MODE" == "pathway" ]]; then
    echo ""
    echo "[2/3] Running Pathway Analysis (using conda env)..."

    # Helper to run R in environment
    RUN_R="${MICROMAMBA} run -p ${PATHWAY_ENV} Rscript"

    cd "${DOWNSTREAM_DIR}/pathway_analysis"

    # 1. Prepare Data
    echo "  -> Running prepare_data.R..."
    ${RUN_R} scripts/prepare_data.R

    # 2. Run GSEA & ORA
    echo "  -> Running Gene Set Enrichment (GSEA)..."
    ${RUN_R} scripts/run_gsea.R
    echo "  -> Running Over-Representation Analysis (ORA)..."
    ${RUN_R} scripts/run_ora.R

    # 3. Running Refactored Pipeline Components
    echo "  -> Running ssGSEA and Visualization..."
    ${RUN_R} scripts/run_refactored_ssgsea.R
    ${RUN_R} scripts/generate_species_contrast_plots.R
    ${RUN_R} scripts/generate_creative_plots.R
    ${RUN_R} scripts/generate_summary_plots.R

    echo "[OK] Pathway Analysis Complete."
else
    echo "[SKIP] Pathway Analysis"
fi

# ------------------------------------------------------------------
# 3. Patient Analysis Comparisons
# ------------------------------------------------------------------
if [[ "$MODE" == "all" || "$MODE" == "patient" ]]; then
    echo ""
    echo "[3/3] Running Patient RNA-seq Downstream Comparisons..."
    
    # Run Cross-Dataset Comparison (NAS/Fibrosis consistency)
    cd "${PROJECT_ROOT}/RNA-seq/patient_RNAseq"
    python3 analysis/downstream/cross_dataset_compare.py

    echo "[OK] Patient Analysis Complete."
else
    echo "[SKIP] Patient Analysis"
fi

echo ""
echo "=========================================================="
echo "   Pipeline Execution Finished"
echo "=========================================================="
