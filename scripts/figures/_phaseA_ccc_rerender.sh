#!/usr/bin/env bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/phaseA_ccc_rerender_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/phaseA_ccc_rerender_%j.out

# Phase A CCC re-render — mega-review A6/A7 remediation.
#
# The LIANA differential producer (fig2_advanced_analyses.py) has been re-run and
# sign-VERIFIED: Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv
# now uses the CORRECTED convention (mean score_diff = -0.303; score_diff>0 == MASLD-enriched).
#
# This script (1) re-runs the CCC consumers that read that CSV (and its downstream
# LR-differential / consensus / edge-atlas chain), then (2) re-renders the
# LIANA-dependent CCC figures.
#
# DEFERRED to Phase B (need the spatial COMMOT/squidpy re-runs):
#   scripts/figures/figS_spatial_ccc_consensus.R, Analysis/Spatial/scripts/30_spatial_ccc_consensus.R
#
# NOTE: do NOT use `set -u` — micromamba's binutils activation hook references an
# unbound ADDR2LINE and trips it (matches run_pub_figures.sh / rerender_fig4_spatial_dependent.sh).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$PROJ"
cd "$PROJ"
mkdir -p "$PROJ/scripts/figures/logs"

# use the env's interpreters directly to avoid PATH ambiguity
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PYBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

echo "=================================================================="
echo "host=$(hostname) start=$(date)"
echo "MASLD_PROJECT_ROOT=$MASLD_PROJECT_ROOT"
echo "liana CSV mtime: $(stat -c '%y' "$PROJ/Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv")"
echo "=================================================================="

# Per-step isolation: one failure must NOT abort the rest.
run_R () {
  local f="$1"
  echo ""
  echo "------------------------------------------------------------------"
  echo "=== [Rscript] $f  ($(date +%H:%M:%S)) ==="
  echo "------------------------------------------------------------------"
  set +e
  "$RBIN" "$f"
  local rc=$?
  set -e
  if [ "$rc" -ne 0 ]; then echo "WARN: $f exited rc=$rc"; else echo "OK: $f"; fi
}

run_PY () {
  local f="$1"
  echo ""
  echo "------------------------------------------------------------------"
  echo "=== [python] $f  ($(date +%H:%M:%S)) ==="
  echo "------------------------------------------------------------------"
  set +e
  "$PYBIN" "$f"
  local rc=$?
  set -e
  if [ "$rc" -ne 0 ]; then echo "WARN: $f exited rc=$rc"; else echo "OK: $f"; fi
}

# ============================================================
# (1) CCC CONSUMERS — recompute on the corrected liana CSV
#     (LR-differential -> consensus -> edge-atlas chain ordered so producers
#      run before downstream consumers: 284 -> 288 -> 253 -> 290)
# ============================================================
run_R "$PROJ/Analysis/SingleCell/scripts/319_nichenet_bulk_reverse_validation.R"
run_R "$PROJ/Analysis/SingleCell/scripts/323_kupffer_to_lam_trajectory.R"
run_R "$PROJ/Analysis/SingleCell/scripts/324_hep_stromal_circuits.R"
run_R "$PROJ/Analysis/SingleCell/scripts/326_hep_metabolic_macrophage_ccc.R"
run_R "$PROJ/Analysis/SingleCell/scripts/329_chromatin_ccc_coupling.R"
run_R "$PROJ/Analysis/SingleCell/scripts/337_sex_stratified_ccc.R"
run_R "$PROJ/Analysis/SingleCell/scripts/339_crossspecies_ccc_conservation.R"
run_R "$PROJ/RNA-seq/85_secretome_plasma_chain.R"
run_R "$PROJ/RNA-seq/284_lr_differential.R"
run_R "$PROJ/RNA-seq/288_lr_consensus_edges.R"
run_R "$PROJ/RNA-seq/253_regulon_lr_cerna_edges.R"
run_PY "$PROJ/RNA-seq/290_edge_annotation_atlas.py"
run_R "$PROJ/Analysis/SingleCell/scripts/compute_translational_priority.R"

# ============================================================
# (2) LIANA-DEPENDENT CCC FIGURES — render after consumers refresh inputs
#     (figS_spatial_ccc_consensus + 30_spatial_ccc_consensus DEFERRED to Phase B)
# ============================================================
cd "$PROJ/scripts/figures"
run_R "$PROJ/scripts/figures/figS_kupffer_lam_trajectory.R"
run_R "$PROJ/scripts/figures/figS_hep_stromal_circuits.R"
run_R "$PROJ/scripts/figures/figS_hep_metabolic_mac.R"
run_R "$PROJ/scripts/figures/figS_chromatin_ccc_coupling.R"
run_R "$PROJ/scripts/figures/figS_crossspecies_ccc.R"
run_R "$PROJ/scripts/figures/figS_sex_stratified_ccc.R"
run_R "$PROJ/scripts/figures/figS_liana_bulk_reverse.R"
run_R "$PROJ/scripts/figures/figS_celltype_biology_overview.R"
run_R "$PROJ/scripts/figures/figS_secretome_plasma_chain.R"
run_R "$PROJ/scripts/figures/regen_liana_crosstalk.R"

echo ""
echo "=================================================================="
echo "=== DONE $(date) ==="
echo "=================================================================="
