#!/bin/bash
# run_gwas_rna_figures.sh — Generate GWAS-RNA-scRNA co-analysis figure panels
#
# Run after Phase A scripts (200-205) have completed.
# Generates panels for Fig 4 and Fig 6.
#
# Usage: bash scripts/figures/run_gwas_rna_figures.sh

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
LOGDIR="${BASE}/RNA-seq/logs/gwas_rna_integration"
FIGDIR="${BASE}/scripts/figures"
mkdir -p "$LOGDIR"

export MASLD_PROJECT_ROOT="$BASE"

echo "=== GWAS-RNA-scRNA Figure Generation ==="
echo "Start: $(date)"

# Fig 4 panels
sbatch --parsable \
  --partition=cpu --cpus-per-task=4 --mem=32G --time=48:00:00 \
  --job-name=fig4_gwas_volcano \
  --output="${LOGDIR}/fig_207_%j.out" \
  --error="${LOGDIR}/fig_207_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${FIGDIR}/fig4_gwas_volcano.R"

# RETIRED 2026-07-07: fig4_expression_pip_scatter.R + fig4_celltype_heritability.R dropped
# from this orchestrator — their outputs (panel_expression_*/panel_celltype_*/panel_geneset_*)
# are stale panels not in the Fig2/FigS2 set; the scripts themselves now early-quit.

# Fig 6 panels — archived 2026-04-22 (6-fig → 5-fig restructure; convergence content absorbed into fig5_convergence)
# Legacy fig6_convergence_sankey.R lives in archive/fig6_ditched_2026-04-22/.

echo ""
echo "=== All figure jobs submitted ==="
echo "Monitor: squeue -u \$USER | grep fig"
