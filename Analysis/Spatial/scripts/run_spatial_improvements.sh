#!/usr/bin/env bash
##############################################################################
# run_spatial_improvements.sh — Submit spatial analysis improvement pipeline
#
# Execution order:
#   Steps 2,3,4 in PARALLEL (03c, 04b, 05b)
#     ↓ --dependency=afterok
#   Step 5: 06_integration.py
#     ↓ --dependency=afterok
#   Step 6: 07_spatial_figures.py
##############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== Spatial Analysis Improvements Pipeline ==="
echo "  Script dir: $SCRIPT_DIR"
echo ""

# Phase 1: Three parallel jobs (03c, 04b, 05b)
JOB_03C=$(sbatch --parsable "$SCRIPT_DIR/03c_validate_deconv.sh")
echo "  Submitted 03c_validate_deconv: SLURM $JOB_03C"

JOB_04B=$(sbatch --parsable "$SCRIPT_DIR/04b_deg_zonation_mapping.sh")
echo "  Submitted 04b_deg_zonation_mapping: SLURM $JOB_04B"

JOB_05B=$(sbatch --parsable "$SCRIPT_DIR/05b_ligand_receptor.sh")
echo "  Submitted 05b_ligand_receptor: SLURM $JOB_05B"

# Phase 2: Integration (depends on all three)
JOB_06=$(sbatch --parsable --dependency=afterok:${JOB_03C}:${JOB_04B}:${JOB_05B} "$SCRIPT_DIR/06_integration.sh")
echo "  Submitted 06_integration: SLURM $JOB_06 (depends on $JOB_03C,$JOB_04B,$JOB_05B)"

# Phase 3: Figures (depends on integration)
JOB_07=$(sbatch --parsable --dependency=afterok:${JOB_06} "$SCRIPT_DIR/07_spatial_figures.sh")
echo "  Submitted 07_spatial_figures: SLURM $JOB_07 (depends on $JOB_06)"

echo ""
echo "=== Pipeline submitted ==="
echo "  Monitor: squeue -u \$USER"
echo "  Logs: $SCRIPT_DIR/../logs/"
echo ""
echo "  Job chain:"
echo "    Phase 1 (parallel): $JOB_03C, $JOB_04B, $JOB_05B"
echo "    Phase 2 (integration): $JOB_06"
echo "    Phase 3 (figures): $JOB_07"
