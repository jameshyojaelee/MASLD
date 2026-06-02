#!/bin/bash
# run_coloc_coworker2.sh — SuSiE-COLOC jobs for coworker2
# Runs 3 GWAS, excluding chr already done or actively running on jameslee's account.
#
# Usage:
#   cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
#   bash src/run_coloc_coworker2.sh
#
# Monitor:
#   squeue -u $(whoami) -n liver-coloc

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

echo "SuSiE-COLOC coworker2 batch — $(date)"
echo ""

# --- Ghouse_HCC: skip chr 4 (done), 6,7,20 (running) ---
# Remaining: 1-3,5,8-19,21-22
SAFE="3,9-17,19,21-22"
RISKY="1-2,5,8,18"
JID1=$(GWAS_NAME=Ghouse_HCC sbatch --parsable --mem=96G  --array=${SAFE}  src/06_susie_coloc_coworker2.sh)
JID2=$(GWAS_NAME=Ghouse_HCC sbatch --parsable --mem=128G --array=${RISKY} src/06_susie_coloc_coworker2.sh)
echo "  Ghouse_HCC: safe=${JID1} risky=${JID2}  (skip chr 4,6,7,20)"

# --- FinnGen_HCC: skip chr 6 (done), 4,7,20 (running) ---
# Remaining: 1-3,5,8-19,21-22
SAFE="3,9-17,19,21-22"
RISKY="1-2,5,8,18"
JID1=$(GWAS_NAME=FinnGen_HCC sbatch --parsable --mem=96G  --array=${SAFE}  src/06_susie_coloc_coworker2.sh)
JID2=$(GWAS_NAME=FinnGen_HCC sbatch --parsable --mem=128G --array=${RISKY} src/06_susie_coloc_coworker2.sh)
echo "  FinnGen_HCC: safe=${JID1} risky=${JID2}  (skip chr 4,6,7,20)"

# --- 2023_36280732_NAFLD_UKBB_EUR: skip chr 20 (done), 4,6,7 (running) ---
# Remaining: 1-3,5,8-19,21-22
SAFE="3,9-17,19,21-22"
RISKY="1-2,5,8,18"
JID1=$(GWAS_NAME=2023_36280732_NAFLD_UKBB_EUR sbatch --parsable --mem=96G  --array=${SAFE}  src/06_susie_coloc_coworker2.sh)
JID2=$(GWAS_NAME=2023_36280732_NAFLD_UKBB_EUR sbatch --parsable --mem=128G --array=${RISKY} src/06_susie_coloc_coworker2.sh)
echo "  2023_36280732_NAFLD_UKBB_EUR: safe=${JID1} risky=${JID2}  (skip chr 4,6,7,20)"

echo ""
echo "Submitted 3 GWAS × 18 chr = 54 tasks (9 chr handled by jameslee's running jobs)"
echo "Monitor: squeue -u \$(whoami) -n liver-coloc"
