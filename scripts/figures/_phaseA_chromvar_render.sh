#!/usr/bin/env bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/phaseA_chromvar_render_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/phaseA_chromvar_render_%j.out

# Phase A chromVAR re-render — mega-review A6 pseudoreplication remediation.
#
# The three chromVAR figure scripts have been repointed (A6.1 + M2) off the
# retired per-CELL Mann-Whitney table (chromvar_tf_activity.csv; 4,832 "sig" hits
# were n-of-cells inflation, treating each cell as an independent replicate of an
# <=18-donor contrast) onto the DONOR-LEVEL limma table:
#   Analysis/ATAC/Human_Multiome/results/chromvar_v2/chromvar_limma_per_ct.csv
# Verified on disk (2026-06-21): repoint live in
#   fig3_chromvar_tfct_dotplot.R:29  -> chromvar_limma_per_ct.csv
#   figS05_gwas_atac.R:186           -> chromvar_limma_per_ct.csv (per-cell = fallback only)
#   fig3_epigenomic_panels.R:241-247 -> donor-level schema shim
# Donor-level limma = 110 sig (adj.P.Val < 0.05) genome-wide, ALL in the
# Low_confidence cell type; 0 among the displayed cell types (incl. Hepatocytes).
#
# HONEST RESULT: these panels render near-NULL (no chromVAR-significant motif in
# any displayed cell type). That is the correct, pseudoreplication-free answer —
# the user may choose to CUT these panels rather than show an all-grey/empty plot.
#
# NOTE: do NOT use `set -u` — micromamba's binutils activation hook references an
# unbound ADDR2LINE and trips it (matches _phaseA_ccc_rerender.sh / run_pub_figures.sh).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$PROJ"
cd "$PROJ"
mkdir -p "$PROJ/scripts/figures/logs"

# use the env's interpreter directly to avoid PATH ambiguity
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

CHROMVAR_CSV="$PROJ/Analysis/ATAC/Human_Multiome/results/chromvar_v2/chromvar_limma_per_ct.csv"

echo "=================================================================="
echo "host=$(hostname) start=$(date)"
echo "MASLD_PROJECT_ROOT=$MASLD_PROJECT_ROOT"
echo "chromVAR donor CSV: $CHROMVAR_CSV"
echo "  mtime: $(stat -c '%y' "$CHROMVAR_CSV" 2>/dev/null || echo MISSING)"
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

# ============================================================
# chromVAR figures — re-render on the donor-level limma table
# ============================================================
cd "$PROJ/scripts/figures"
run_R "$PROJ/scripts/figures/fig3_chromvar_tfct_dotplot.R"
run_R "$PROJ/scripts/figures/fig3_epigenomic_panels.R"
run_R "$PROJ/scripts/figures/figS05_gwas_atac.R"

echo ""
echo "=================================================================="
echo "=== DONE $(date) ==="
echo "=================================================================="
