#!/bin/bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Cas13_Library_Design/logs/figs_v8_%j.out
#SBATCH --error=Cas13_Library_Design/logs/figs_v8_%j.err

# Re-render ALL Cas13 library figures for v9.
#   - human + mouse DE now both use the limma-voom engine
#       human  = canonical_deg_results.csv (limma-voom QW C2)
#       mouse  = per_diet_cas13/*_de_results.csv (limma-voom + ashr)  [unchanged]
#   - no miRNA anywhere (Cas13 cannot knock down miRNA)
#   - subtitles removed, titles simplified
# Runs each script independently; one failure does NOT block the rest.
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Fig 10 (lfc-vs-tpm strategy) and Fig 12c (tpm recovery) retired 2026-06-29 (TPM gate
# replaced by sc-hepatocyte CPM); their scripts moved to *_legacy/. Expression-cutoff
# sensitivity now lives solely in fig13_cpm_cutoff_sensitivity.R (13a-b).
SCRIPTS=(
  "scripts/figures/figS_cas13_library_diet.R"                    # 01, 02
  "scripts/figures/figS_cas13_library_ortholog_lncrna.R"        # 04  (miRNA removed)
  "scripts/figures/figS_cas13_library_human.R"                  # 05
  "Cas13_Library_Design/scripts/benchmark_cutoffs.R"            # 07a-d
  "scripts/figures/figS_raw_vs_shrunk_patient_concordance.R"    # 08a-d
  "scripts/figures/figS_cas13_library_options.R"               # 11a-c
  "scripts/figures/fig12_positive_control_recovery.R"          # 12b
  "scripts/figures/figS_cas13_library_core_vs_mash_overlap.R"  # core_vs_cohort_euler / evidence_upset / tier_composition / mouse_hep_substrate
  "scripts/figures/fig13_cpm_cutoff_sensitivity.R"             # 13a-b
)

declare -a RESULT
for s in "${SCRIPTS[@]}"; do
  echo "==================== RENDER: $s ===================="
  if Rscript "$s"; then
    RESULT+=("PASS  $s")
  else
    RESULT+=("FAIL  $s (exit $?)")
  fi
done

echo ""
echo "==================== SUMMARY ===================="
printf '%s\n' "${RESULT[@]}"
