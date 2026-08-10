#!/usr/bin/env bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_evidence_classes_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_evidence_classes_%j.out

# Re-run scripts/figures/fig1_fig4_evidence_classes.R UNCHANGED against the frozen
# 2026-07-15-r2 manuscript release, regenerating the three evidence-class Figure 4
# panels (fig4a positive rates / fig4b adjusted OR / fig4c convergence comparison)
# into FIG4_DIR root.
#
# The script default is already MANUSCRIPT_RELEASE_ID=2026-07-15-r2 (line 16); the
# export below pins it explicitly so an inherited environment cannot silently point
# the run at the stale 2026-07-10-r1 release.
#
# COLLISION NOTE: this script ALSO writes two Figure 1 panels
#   figures/main/fig1_atlas_overview/fig1b_complementary_maps.pdf   (line 86-87)
#   figures/main/fig1_atlas_overview/fig1c_genetic_trait_scope.pdf  (line 112-113)
# Both are backed up before submission and their hashes are printed below so the
# re-run can be confirmed non-destructive.

set -eo pipefail
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MANUSCRIPT_RELEASE_ID=2026-07-15-r2

echo "host=$(hostname) start=$(date) cwd=$(pwd)"
echo "release=$MANUSCRIPT_RELEASE_ID"

echo "=== [Rscript] fig1_fig4_evidence_classes.R ($(date +%H:%M:%S)) ==="
set +e
"$RBIN" scripts/figures/fig1_fig4_evidence_classes.R
rc=$?
set -e
echo "rc=$rc"

echo "=== Figure 4 evidence-class outputs ==="
ls -la --time-style=full-iso figures/main/fig4_validation/*.pdf || true
sha256sum figures/main/fig4_validation/*.pdf || true

echo "=== Figure 1 collision check (expect unchanged content) ==="
sha256sum figures/main/fig1_atlas_overview/fig1b_complementary_maps.pdf \
          figures/main/fig1_atlas_overview/fig1c_genetic_trait_scope.pdf || true

echo "=== archived r2 reference hashes (for diff) ==="
sha256sum figures/main/fig4_validation/_legacy/2026-08-06_preconsolidation/fig4a_evidence_class_positive_rates.pdf \
          figures/main/fig4_validation/_legacy/2026-08-06_preconsolidation/fig4b_evidence_class_adjusted_or.pdf \
          figures/main/fig4_validation/_legacy/2026-08-06_preconsolidation/fig4c_convergence_comparison.pdf || true

echo "=== DONE $(date) ==="
