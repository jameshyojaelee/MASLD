#!/usr/bin/env bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/phaseA_ccc_figs_only_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/phaseA_ccc_figs_only_%j.out

# Phase A CCC figure re-render (figures ONLY — the consumers already re-ran
# successfully under job 17513236; secretome_chain.csv regenerated 2026-06-21 11:00).
# Fixes the cwd bug in _phaseA_ccc_rerender.sh:94 (it cd'd into scripts/figures, but the
# figS_*.R scripts source the PROJECT-ROOT-relative "scripts/figures/publication_theme.R").
# Therefore we render from the PROJECT ROOT.
# DEFERRED to Phase B (need spatial COMMOT/squidpy re-runs): figS_spatial_ccc_consensus.R.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$PROJ"
cd "$PROJ"                      # <-- project root, NOT scripts/figures
mkdir -p "$PROJ/scripts/figures/logs"
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

echo "host=$(hostname) start=$(date)  cwd=$(pwd)"
echo "liana CSV mtime: $(stat -c '%y' "$PROJ/Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv")"

run_R () {
  echo ""; echo "=== [Rscript] $1  ($(date +%H:%M:%S)) ==="
  set +e; "$RBIN" "$1"; local rc=$?; set -e
  if [ "$rc" -ne 0 ]; then echo "WARN: $1 exited rc=$rc"; else echo "OK: $1"; fi
}

for f in \
  figS_kupffer_lam_trajectory figS_hep_stromal_circuits figS_hep_metabolic_mac \
  figS_chromatin_ccc_coupling figS_crossspecies_ccc figS_sex_stratified_ccc \
  figS_liana_bulk_reverse figS_celltype_biology_overview figS_secretome_plasma_chain \
  regen_liana_crosstalk ; do
  run_R "$PROJ/scripts/figures/${f}.R"
done

echo ""; echo "=== DONE $(date) ==="
