#!/bin/bash
# run_coloc_coworker3.sh — SuSiE-COLOC jobs for coworker3
# Runs 3 GWAS, excluding chr already done.
#
# Usage:
#   cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
#   bash src/run_coloc_coworker3.sh
#
# Monitor:
#   squeue -u $(whoami) -n liver-GWAS

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

echo "SuSiE-COLOC coworker3 batch — $(date)"
echo ""

# --- PanUKBB_CSA_ALT: skip chr 4,6,7,20 (done) ---
# Remaining: 1-3,5,8-19,21-22
SAFE="3,9-17,19,21-22"
RISKY="1-2,5,8,18"
JID1=$(GWAS_NAME=PanUKBB_CSA_ALT sbatch --parsable --mem=96G  --array=${SAFE}  src/06_susie_coloc_coworker3.sh)
JID2=$(GWAS_NAME=PanUKBB_CSA_ALT sbatch --parsable --mem=128G --array=${RISKY} src/06_susie_coloc_coworker3.sh)
echo "  PanUKBB_CSA_ALT: safe=${JID1} risky=${JID2}  (skip chr 4,6,7,20)"

# --- PanUKBB_CSA_AST: skip chr 4,6,7,20 (done) ---
# Remaining: 1-3,5,8-19,21-22
SAFE="3,9-17,19,21-22"
RISKY="1-2,5,8,18"
JID1=$(GWAS_NAME=PanUKBB_CSA_AST sbatch --parsable --mem=96G  --array=${SAFE}  src/06_susie_coloc_coworker3.sh)
JID2=$(GWAS_NAME=PanUKBB_CSA_AST sbatch --parsable --mem=128G --array=${RISKY} src/06_susie_coloc_coworker3.sh)
echo "  PanUKBB_CSA_AST: safe=${JID1} risky=${JID2}  (skip chr 4,6,7,20)"

# --- 2019_31311600_NAFLD_EUR: skip chr 13,18,21,22 (done) ---
# Remaining: 1-12,14-17,19-20
SAFE="3-4,6-7,9-12,14-17,19-20"
RISKY="1-2,5,8"
JID1=$(GWAS_NAME=2019_31311600_NAFLD_EUR sbatch --parsable --mem=96G  --array=${SAFE}  src/06_susie_coloc_coworker3.sh)
JID2=$(GWAS_NAME=2019_31311600_NAFLD_EUR sbatch --parsable --mem=128G --array=${RISKY} src/06_susie_coloc_coworker3.sh)
echo "  2019_31311600_NAFLD_EUR: safe=${JID1} risky=${JID2}  (skip chr 13,18,21,22)"

echo ""
echo "Submitted 3 GWAS × 18 chr = 54 tasks"
echo "Monitor: squeue -u \$(whoami) -n liver-GWAS"
