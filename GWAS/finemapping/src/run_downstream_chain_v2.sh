#!/bin/bash
# run_downstream_chain_v2.sh
# Waits for chr1 eQTL SuSiE, submits chr1 COLOC at 128G,
# then chains downstream: Script 07 → 27a → 75 → 48/217
# All using afterok on EVERY COLOC job currently in the queue.

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

GWAS_LIST=(
  2019_31311600_NAFLD_EUR 2020_32298765_NAFLD_EUR 2021_34128465_PDFF_EUR
  2021_34841290_NAFLD_EUR 2021_34957434_PDFF_EUR 2022_36402844_PDFF_EUR
  2023_36280732_NAFLD_deCode_EUR 2023_36280732_NAFLD_Intermountain_EUR
  2023_36280732_NAFLD_UKBB_EUR UKBB_ALT UKBB_AST UKBB_GGT
  FinnGen_NAFLD FinnGen_NASH FinnGen_HCC Ghouse_Cirrhosis Ghouse_HCC
  2020_32514122_Cirrhosis_EAS 2020_32514122_HCC_EAS
  BBJ_ALT BBJ_AST BBJ_GGT
  PanUKBB_AFR_ALT PanUKBB_AFR_AST PanUKBB_AFR_GGT
  PanUKBB_CSA_ALT PanUKBB_CSA_AST PanUKBB_CSA_GGT
)

echo "============================================"
echo "Downstream Chain v2 — $(date)"
echo "============================================"

# ── Step 1: Wait for chr1 eQTL SuSiE, then submit chr1 COLOC at 128G ──
echo "--- Waiting for chr1 eQTL SuSiE ---"
while true; do
  count=$(find results/eqtl_susie/chr1/ -name "*.rds" 2>/dev/null | wc -l)
  workers=$(squeue --me -n eqtl_susie --format="%i" 2>/dev/null | grep "_1$" | wc -l)
  if [ $workers -eq 0 ] || [ $count -ge 2048 ]; then
    echo "  chr1 READY: ${count}/2048 RDS, ${workers} workers"
    break
  fi
  echo "  chr1: ${count}/2048 (${workers} workers) — waiting 2 min..."
  sleep 120
done

echo "Submitting chr1 COLOC at 128G..."
CHR1_JOBS=""
for gwas in "${GWAS_LIST[@]}"; do
  JID=$(GWAS_NAME="${gwas}" sbatch --parsable --mem=128G --array=1 src/06_susie_coloc.sh)
  CHR1_JOBS="${CHR1_JOBS}${CHR1_JOBS:+,}${JID}"
done
echo "  Submitted 28 chr1 COLOC jobs at 128G"

# ── Step 2: Collect ALL COLOC job IDs currently in queue ──
sleep 5
ALL_COLOC=$(squeue --me -n susie_coloc --format="%i" 2>/dev/null | grep -o "^[0-9]*" | sort -un | tr '\n' ':' | sed 's/:$//')
echo ""
echo "All COLOC parent job IDs for dependency: ${ALL_COLOC}"
N_PARENTS=$(echo "$ALL_COLOC" | tr ':' '\n' | wc -l)
echo "  (${N_PARENTS} parent array jobs)"

# ── Step 3: Script 07 — combine all COLOC results ──
echo ""
echo "--- Submitting Script 07 (combine COLOC) ---"
JOB_07=$(sbatch --parsable \
  --dependency=afterok:${ALL_COLOC} \
  --partition=cpu --mem=32G --time=4:00:00 \
  --job-name=combine_coloc \
  --output=logs/combine_coloc_%j.out \
  --error=logs/combine_coloc_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate finemapping && Rscript src/07_combine_susie_coloc.R")
echo "  Job 07: ${JOB_07} (afterok: ALL COLOC)"

# ── Step 4: Script 27a — rebuild atlas ──
echo ""
echo "--- Submitting Script 27a (atlas assembly) ---"
JOB_27a=$(sbatch --parsable \
  --dependency=afterok:${JOB_07} \
  --partition=cpu --mem=64G --time=8:00:00 \
  --job-name=atlas_27a \
  --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/atlas_27a_%j.out \
  --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/atlas_27a_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate rnaseq && Rscript RNA-seq/27a_assemble_evidence_atlas.R")
echo "  Job 27a: ${JOB_27a} (afterok: ${JOB_07})"

# ── Step 5: Script 75 — integrate causal overhaul ──
echo ""
echo "--- Submitting Script 75 (causal overhaul) ---"
JOB_75=$(sbatch --parsable \
  --dependency=afterok:${JOB_27a} \
  --partition=cpu --mem=64G --time=8:00:00 \
  --job-name=causal_75 \
  --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/causal_75_%j.out \
  --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/causal_75_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate rnaseq && Rscript RNA-seq/75_integrate_causal_overhaul.R")
echo "  Job 75: ${JOB_75} (afterok: ${JOB_27a})"

# ── Step 6: Script 48 + Script 217 (parallel after 75) ──
echo ""
echo "--- Submitting Script 48 (cross-ancestry) ---"
JOB_48=$(sbatch --parsable \
  --dependency=afterok:${JOB_75} \
  --partition=cpu --mem=32G --time=4:00:00 \
  --job-name=cross_ancestry_48 \
  --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cross_ancestry_48_%j.out \
  --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cross_ancestry_48_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate rnaseq && Rscript RNA-seq/48_cross_ancestry_replication.R")
echo "  Job 48: ${JOB_48} (afterok: ${JOB_75})"

echo ""
echo "--- Submitting Script 217 (stratified causal) ---"
JOB_217=$(sbatch --parsable \
  --dependency=afterok:${JOB_75} \
  --partition=cpu --mem=32G --time=4:00:00 \
  --job-name=stratified_217 \
  --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/stratified_217_%j.out \
  --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/stratified_217_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate rnaseq && Rscript RNA-seq/217_stratified_causal_atlas.R")
echo "  Job 217: ${JOB_217} (afterok: ${JOB_75})"

echo ""
echo "============================================"
echo "FULL DEPENDENCY CHAIN:"
echo "  All COLOC (chr1@128G + chr2@128G + chr3@96G + chr4-22 split)"
echo "    └→ Script 07  : ${JOB_07}"
echo "        └→ Script 27a : ${JOB_27a}"
echo "            └→ Script 75  : ${JOB_75}"
echo "                ├→ Script 48  : ${JOB_48}"
echo "                └→ Script 217 : ${JOB_217}"
echo "============================================"
