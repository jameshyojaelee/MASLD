#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/regen_fig2_fig5_panels_%j.out
#SBATCH --error=scripts/figures/logs/regen_fig2_fig5_panels_%j.err
# LEGACY PRE-RESOURCE regeneration path. It writes retired Figure 5 convergence
# panels and therefore fails closed by default. The directory name
# `fig5_convergence` is retained only for backward compatibility; current
# Figure 5 is the MASLD Gene Catalog.
# Excludes: fig2_progression.R, assemble_*, fig5_causal_architecture.R,
#   fig5_translation.R, fig5_convergence_legacy.R, *_compact. Composite byproducts
#   are purged at the end (explicit filenames only, to protect real panels).
if [[ "${ALLOW_LEGACY_FIGURE_REGEN:-false}" != "true" ]]; then
  echo "REFUSED: this script regenerates retired convergence panels."
  echo "Use Plan 60 for the Resource release; set ALLOW_LEGACY_FIGURE_REGEN=true only for provenance-only regeneration."
  exit 64
fi

set +e
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE" || exit 1
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
SF=scripts/figures

run() { echo "===== $1 ====="; Rscript "$SF/$1" 2>&1 | tail -3; echo "[exit $1=${PIPESTATUS[0]}]"; }

echo "########## FIG 2 PANELS (standalone; fig2_progression.R composite EXCLUDED) ##########"
for s in bayesprism_transitions canonical_vs_perstudy_upset cascade_degs \
         celltype_bulk_attribution celltype_cascade celltype_concordance \
         fib_stage_degs fib_stage_upset \
         fib_stage_vs_ctrl_degs fib_stage_vs_ctrl_upset hkdc1_module24 \
         nas_fib_grid nas_stage_degs nas_stage_upset \
         nas_stage_vs_ctrl_degs nas_stage_vs_ctrl_upset network_communities; do
  [ -f "$SF/$s.R" ] && run "$s.R"
done
for extra in ccc_v3_panels.R fig2_chromatin_cascade.R hotspot_module_geneset_correlation.R gen_scrna_umap_embeddable.R; do
  [ -f "$SF/$extra" ] && run "$extra"
done

echo "########## FIG 5 PANELS (assemble/causal_arch/translation/legacy EXCLUDED) ##########"
run fig5_convergence.R                              # -> panels/fig6_therapeutic_axes.pdf (canonical)
run fig5_convergence_v3.R                           # -> panels/fig5b.pdf
# panel 5e (TF convergence scatter + 4-way survival lollipop) CUT 2026-06-19:
#   banned lollipop + oversold refuted 4-way claim + null SCENIC+ axis.
#   Source scripts moved to scripts/figures/_legacy/.

echo "########## PURGE composite byproducts (explicit names only) ##########"
for c in fig2_progression.pdf fig2_progression_sex.pdf fig5_convergence.pdf fig5_composite.pdf fig5_composite_v2.pdf fig1_compact.pdf; do
  find figures/main -name "$c" -delete -print 2>/dev/null
done

echo "########## RESULTING PANELS in figures/main ##########"
echo "-- fig3_RNAseq --"; find figures/main/fig3_RNAseq -name "*.pdf" 2>/dev/null | sed 's#.*/##' | sort
echo "-- fig5_convergence --";     find figures/main/fig5_convergence     -name "*.pdf" 2>/dev/null | sed 's#.*/##' | sort
echo "DONE"
