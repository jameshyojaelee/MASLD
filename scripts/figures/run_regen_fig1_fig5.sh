#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/regen_fig1_fig5_%j.out
#SBATCH --error=scripts/figures/logs/regen_fig1_fig5_%j.err
# Regenerate Fig 1 (atlas overview) + Fig 5 (convergence + drug targets) against
# the rebuilt 389-col atlas / refreshed convergence_evidence.csv. RUN-ONLY — does
# not modify any figure script. set +e so one panel failing doesn't abort the rest.
set +e
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE" || exit 1
mkdir -p scripts/figures/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

step() { echo "=================== $* ==================="; }

step "FIG 1 — py panels"
python scripts/figures/fig_sunburst.py;                 echo "[exit fig_sunburst=$?]"
python scripts/figures/fig_convergence_wheel.py --no-title; echo "[exit convergence_wheel=$?]"
python scripts/figures/fig_bulkrna_matrix.py --panel;   echo "[exit bulkrna_matrix=$?]"
step "FIG 1 — compact assembler"
Rscript scripts/figures/fig1_compact.R;                 echo "[exit fig1_compact=$?]"

step "FIG 3 — canonical DEG volcano (limma-voom C2; moved from fig1)"
Rscript scripts/figures/fig3_deg_volcano.R;             echo "[exit fig3_deg_volcano=$?]"

step "FIG 1 — RNA-seq concordance data (05_concordance)"
INT=RNA-seq/Human/Patient_Cohorts/analysis/integration
Rscript "${INT}/scripts/mega_validation/05_concordance.R"; echo "[exit 05_concordance=$?]"

step "FIG S — mega-validation alluvial + chord (figS_mega_validation)"
Rscript scripts/figures/figS_mega_validation.R;         echo "[exit figS_mega_validation=$?]"

step "FIG 5 — panel 5a (multi-evidence matrix)"
Rscript scripts/figures/fig5_convergence.R;             echo "[exit fig5_convergence(5a)=$?]"
step "FIG 5 — panel 5b (sources-active + recovery)"
Rscript scripts/figures/fig5_convergence_v3.R;          echo "[exit fig5_convergence_v3(5b)=$?]"
# Composite assembly removed 2026-06-20 (PI directive): individual panels only.

step "OUTPUT PDFs (mtime)"
ls -l --time-style=+%H:%M \
      figures/main/fig1_atlas_overview/fig1_compact.pdf \
      figures/main/fig3_RNAseq/panels/deg_volcano.pdf \
      figures/supplementary/figS_methods_validation/mega_validation/panels/panelF_method_overlap_alluvial.pdf \
      figures/supplementary/figS_methods_validation/mega_validation/panels/panelF_method_overlap_chord.pdf \
      figures/main/fig5_convergence/panels/fig5{a,b}.pdf 2>/dev/null
echo "DONE regen"
