#!/usr/bin/env bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_AB_rerender_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_AB_rerender_%j.out

# Re-render Fig4 Panel B after the 56 motif-disruption producer refresh
# (refreshed motif_in_disease_regulon flag against the current 48-TF set + dedup).
# Panel B reads motif_disruption_scores.csv. Run from PROJECT ROOT (the figS scripts
# source the project-root-relative scripts/figures/publication_theme.R).
# NOTE: Panel A (disease_master_regulators.R) retired 2026-07-07 — archived to
# scripts/figures/archive/; do not re-add to the loop below.
set -eo pipefail
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "host=$(hostname) start=$(date) cwd=$(pwd)"
echo "motif CSV mtime: $(stat -c '%y' GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv)"

for f in fig4_variant_motif_disruption; do
  echo "=== [Rscript] $f  ($(date +%H:%M:%S)) ==="
  set +e; "$RBIN" "scripts/figures/${f}.R"; rc=$?; set -e
  if [ "$rc" -ne 0 ]; then echo "WARN: $f exited rc=$rc"; else echo "OK: $f"; fi
done
echo "=== DONE $(date) ==="
