#!/bin/bash
# run_stratified_causal_full.sh — Re-run full stratified causal pipeline (207-217)
# Modules A-D run in parallel; figures chain after their analysis scripts; 217 runs last.
# Usage: bash RNA-seq/run_stratified_causal_full.sh

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}/RNA-seq"
mkdir -p logs

echo "=== Stratified Causal Pipeline (207-217) ==="
echo "Start: $(date)"

# ── Module A: Sex-Stratified COLOC (207→208→209) ──
JOB207=$(sbatch --parsable run_207_sex_stratified_coloc.sh)
JOB208=$(sbatch --parsable run_208_subtype_coloc.sh)
JOB209=$(sbatch --parsable --dependency=afterok:${JOB207}:${JOB208} run_209_sex_subtype_figures.sh)
echo "  Module A: 207=${JOB207}, 208=${JOB208}, 209=${JOB209}"

# ── Module B: Progression-Stratified COLOC (210→210b,211→212) ──
JOB210=$(sbatch --parsable run_210_progression_coloc.sbatch)
JOB210B=$(sbatch --parsable --dependency=afterok:${JOB210} run_210b.sbatch)
JOB211=$(sbatch --parsable --dependency=afterok:${JOB210} run_211_progression_driver_genetics.sbatch)
JOB212=$(sbatch --parsable --dependency=afterok:${JOB210}:${JOB211} run_212_progression_figures.sbatch)
echo "  Module B: 210=${JOB210}, 210b=${JOB210B}, 211=${JOB211}, 212=${JOB212}"

# ── Module C: Spatial-Stratified COLOC (213→214) ──
JOB213=$(sbatch --parsable run_213_spatial_coloc.sbatch)
JOB214=$(sbatch --parsable --dependency=afterok:${JOB213} run_214_spatial_figures.sbatch)
echo "  Module C: 213=${JOB213}, 214=${JOB214}"

# ── Module D: Pharmacogenomic Mapping (215→216) ──
JOB215=$(sbatch --parsable run_215_pharmacogenomic.sh)
JOB216=$(sbatch --parsable --dependency=afterok:${JOB215} run_216_pharma_figures.sbatch)
echo "  Module D: 215=${JOB215}, 216=${JOB216}"

# ── Module Integration (217: depends on all analysis scripts) ──
JOB217=$(sbatch --parsable \
  --dependency=afterok:${JOB207}:${JOB208}:${JOB210}:${JOB211}:${JOB213}:${JOB215} \
  run_217_stratified_atlas.sbatch)
echo "  Integration: 217=${JOB217}"

echo ""
echo "  Total: 12 jobs submitted (4 modules parallel + figures + integration)"
echo "  Monitor: squeue -u \$(whoami) | grep -E '20[7-9]|21[0-7]|strat'"
echo "  Expected wall-clock: ~5-6 hours"
echo "=== $(date) ==="
