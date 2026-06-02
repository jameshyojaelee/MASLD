#!/bin/bash
#SBATCH --job-name=concordance_2c_5
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/phase2c_5_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/phase2c_5_%j.err

set -euo pipefail

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT_DIR="${MASLD_PROJECT_ROOT}/Analysis/Cross_Species_Concordance/scripts"

# Track pass/fail for summary
declare -A STATUS

PHASES=(
  "02c_ora_concordance.R"
  "03a_wgcna_module_preservation.R"
  "03b_tf_activity_concordance.R"
  "04_concordance_visualizations.R"
  "05_unified_concordance_atlas.R"
)

for script in "${PHASES[@]}"; do
  echo "============================================"
  echo "Running: ${script}"
  echo "Started: $(date)"
  echo "============================================"
  if Rscript "${SCRIPT_DIR}/${script}"; then
    STATUS["${script}"]="PASS"
    echo ">>> ${script}: PASS ($(date))"
  else
    STATUS["${script}"]="FAIL (exit $?)"
    echo ">>> ${script}: FAIL ($(date))"
  fi
  echo ""
done

echo "============================================"
echo "SUMMARY"
echo "============================================"
for script in "${PHASES[@]}"; do
  printf "  %-40s %s\n" "${script}" "${STATUS[${script}]}"
done
echo ""
echo "Finished: $(date)"
