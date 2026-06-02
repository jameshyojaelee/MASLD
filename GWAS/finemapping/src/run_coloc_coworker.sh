#!/bin/bash
# run_coloc_coworker.sh — SuSiE-COLOC jobs for coworker
# Runs 3 GWAS that are at the back of jameslee's queue.
#
# Usage:
#   cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
#   bash src/run_coloc_coworker.sh
#
# Monitor:
#   squeue -u $(whoami) -n MASLD-james

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

GWAS_LIST=(
  FinnGen_NASH
  2023_36280732_NAFLD_deCode_EUR
  2023_36280732_NAFLD_Intermountain_EUR
)

# Safe chr: 96G. Risky chr (large LD blocks): 128G.
SAFE_CHR="3-4,6-7,9-17,19-22"
RISKY_CHR="1-2,5,8,18"

echo "SuSiE-COLOC coworker batch — $(date)"
echo "GWAS: ${GWAS_LIST[*]}"
echo ""

for gwas in "${GWAS_LIST[@]}"; do
  JID_SAFE=$(GWAS_NAME="${gwas}" sbatch --parsable \
    --mem=96G --array=${SAFE_CHR} src/06_susie_coloc_coworker.sh)
  JID_RISKY=$(GWAS_NAME="${gwas}" sbatch --parsable \
    --mem=128G --array=${RISKY_CHR} src/06_susie_coloc_coworker.sh)
  echo "  ${gwas}: safe=${JID_SAFE} risky=${JID_RISKY}"
done

echo ""
echo "Submitted 3 GWAS × 22 chr = 66 tasks"
echo "Monitor: squeue -u \$(whoami) -n MASLD-james"
