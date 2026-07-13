#!/bin/bash -l
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/coloc_figs_refresh_%j.out
#SBATCH --error=scripts/figures/logs/coloc_figs_refresh_%j.err
# Regenerate the 11 ATLAS-INDEPENDENT genetics/COLOC figures on the 2026-07-05
# MVP 50-GWAS COLOC master (gene_level_coloc.csv / susie_coloc_all_gwas.csv).
# The 7 atlas/convergence-dependent figures (Fig2F, Fig5A, figS09 locus-zooms,
# figS_convergence_evidence, convergence_evidence_matrix) are DEFERRED until the
# concurrent-session atlas resync (C10: 27a -> 75 -> 217) lands.
set -o pipefail   # NOT -u: micromamba/conda activate scripts reference unbound vars (ADDR2LINE)
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPTS=(
  fig2_coloc_counts.R
  fig2_finemap_cascade.R
  fig2_coloc_summary.R
  fig2_finemap_resolution.R
  fig2_finemap_headtohead.R
  fig3b_ancestry_coloc_bars.R
  fig3g_cross_ancestry_pp4.R
  fig3a_hybrid.R
  figS04_coloc.R
  figS09_multi_ancestry_coloc.R
  fig4_celltype_coloc_heatmap.R
)

ok=0; fail=0; failed_list=""
for s in "${SCRIPTS[@]}"; do
  echo "======================== RUN $s ($(date '+%H:%M:%S')) ========================"
  if Rscript "scripts/figures/$s" 2>&1; then
    echo "RESULT $s = OK"
    ok=$((ok+1))
  else
    echo "RESULT $s = FAIL (exit $?)"
    fail=$((fail+1)); failed_list="$failed_list $s"
  fi
done
echo "================================================================"
echo "SUMMARY: $ok OK, $fail FAIL"
[ -n "$failed_list" ] && echo "FAILED:$failed_list"
echo "================================================================"
