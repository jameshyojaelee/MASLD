#!/bin/bash
#SBATCH --job-name=cas13figs
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Cas13_Library_Design/logs/figs_v6_%j.out
#SBATCH --error=Cas13_Library_Design/logs/figs_v6_%j.err

# Re-render all Cas13 library figures for v6 (limma-voom+metafor source).
# Runs each script independently; one failure does NOT block the rest.
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Dependency order: compare_lfc_vs_tpm writes CSVs consumed by figS_..._lfc_vs_tpm.
SCRIPTS=(
  "Cas13_Library_Design/scripts/compare_lfc_vs_tpm_strategy.R"
  "scripts/figures/figS_cas13_library_lfc_vs_tpm.R"
  "scripts/figures/figS_cas13_library_human.R"
  "scripts/figures/figS_cas13_library_tpm_by_diet.R"
  "scripts/figures/figS_cas13_library_tpm_by_diet_ashr04.R"
  "scripts/figures/figS_cas13_library_ortholog_lncrna.R"
  "scripts/figures/figS_cas13_library_diet.R"
  "Cas13_Library_Design/scripts/benchmark_cutoffs.R"
  "scripts/figures/fig12_positive_control_recovery.R"
  "scripts/figures/figS_cas13_library_options.R"
  "scripts/figures/figS_cas13_library_core_vs_mash_overlap.R"
  "scripts/figures/figS_raw_vs_shrunk_patient_concordance.R"
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
