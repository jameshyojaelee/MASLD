#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/regen_fig2_fig5_panels_%j.out
#SBATCH --error=scripts/figures/logs/regen_fig2_fig5_panels_%j.err
# Regenerate INDIVIDUAL PANELS for fig2_progression_sex + fig5_convergence on the
# current -s2 data, into figures/main/. PANELS ONLY — NO composite/assembler runs.
# Excludes: fig2_progression.R, assemble_*, fig5_causal_architecture.R,
#   fig5_translation.R, fig5_convergence_legacy.R, *_compact. Composite byproducts
#   are purged at the end (explicit filenames only, to protect real panels).
set +e
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE" || exit 1
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
SF=scripts/figures

run() { echo "===== $1 ====="; Rscript "$SF/$1" 2>&1 | tail -3; echo "[exit $1=${PIPESTATUS[0]}]"; }

echo "########## FIG 2 PANELS (standalone; fig2_progression.R composite EXCLUDED) ##########"
for s in "$SF"/fig2_panel_*.R; do run "$(basename "$s")"; done
for extra in fig2_ccc_v3_panels.R fig2_chromatin_cascade.R fig2_hotspot_module_geneset_correlation.R gen_fig2e_embeddable.R; do
  [ -f "$SF/$extra" ] && run "$extra"
done

echo "########## FIG 5 PANELS (assemble/causal_arch/translation/legacy EXCLUDED) ##########"
run fig5_convergence.R                              # -> panels/fig5a.pdf (guarded for empty regulons)
run fig5_convergence_v3.R                           # -> panels/fig5b.pdf
run fig5_panel_tf_convergence_atac_rna_coloc.R      # -> fig5_tf_convergence_scatter.pdf
run fig5_panel_D_tf_4way_survival_lollipop.R        # -> fig5_tf_4way_survival_lollipop.pdf

echo "########## PURGE composite byproducts (explicit names only) ##########"
for c in fig2_progression.pdf fig2_progression_sex.pdf fig5_convergence.pdf fig5_composite.pdf fig1_compact.pdf; do
  find figures/main -name "$c" -delete -print 2>/dev/null
done

echo "########## RESULTING PANELS in figures/main ##########"
echo "-- fig2_progression_sex --"; find figures/main/fig2_progression_sex -name "*.pdf" 2>/dev/null | sed 's#.*/##' | sort
echo "-- fig5_convergence --";     find figures/main/fig5_convergence     -name "*.pdf" 2>/dev/null | sed 's#.*/##' | sort
echo "DONE"
