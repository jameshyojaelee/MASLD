#!/bin/bash
# run_gwas_atac.sh
# Orchestrates the GWAS-ATAC variant overlap pipeline
# Chains: 04 (aggregate) → 55 (overlap) → 56 (motif) → 57 (figures/atlas)
# Usage: bash run_gwas_atac.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FM_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="${FM_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "============================================"
echo "GWAS-ATAC Regulatory Variant Pipeline"
echo "Script dir: $SCRIPT_DIR"
echo "Log dir:    $LOG_DIR"
echo "============================================"

# Check if aggregated results exist; if not, run Script 04 first
if [ ! -f "${FM_DIR}/results/combined_finemapping.csv" ]; then
    echo ""
    echo "--- Step 0: Aggregating fine-mapping results (Script 04) ---"
    JOB0=$(sbatch --parsable "${SCRIPT_DIR}/run_04.sbatch")
    echo "  Job 0 (aggregate): $JOB0"
    DEP0="--dependency=afterok:${JOB0}"
else
    echo "combined_finemapping.csv exists — skipping Script 04"
    DEP0=""
fi

# Step 1: Script 55 — Variant-Peak Overlap
echo ""
echo "--- Step 1: Variant-Peak Overlap (Script 55) ---"
JOB1=$(sbatch --parsable --job-name=gwasatac ${DEP0} "${SCRIPT_DIR}/run_55.sbatch")
echo "  Job 1 (overlap): $JOB1"

# Step 2: Script 56 — Motif Disruption
echo ""
echo "--- Step 2: Motif Disruption (Script 56) ---"
JOB2=$(sbatch --parsable --job-name=gwasatac --dependency=afterok:${JOB1} "${SCRIPT_DIR}/run_56.sbatch")
echo "  Job 2 (motif): $JOB2"

# Step 3: Script 57 — Figures & Atlas Integration
echo ""
echo "--- Step 3: Figures & Atlas Integration (Script 57) ---"
JOB3=$(sbatch --parsable --job-name=gwasatac --dependency=afterok:${JOB2} "${SCRIPT_DIR}/run_57.sbatch")
echo "  Job 3 (figures): $JOB3"

echo ""
echo "============================================"
echo "Pipeline submitted:"
if [ -n "$DEP0" ]; then
    echo "  04_aggregate  → $JOB0"
fi
echo "  55_overlap    → $JOB1"
echo "  56_motif      → $JOB2"
echo "  57_figures    → $JOB3"
echo ""
echo "Monitor: squeue -u \$USER -n gwas_atac"
echo "Logs:    $LOG_DIR"
echo "============================================"
