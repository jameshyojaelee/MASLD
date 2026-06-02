#!/bin/bash
# run_downstream_chain.sh
# Submits chr1-3 COLOC + all downstream jobs with strict dependency chain.
# Waits for eQTL SuSiE to finish, submits chr1-3, then chains:
#   All COLOC (chr4-22 already running + chr1-3) → Script 07 → Script 27a → Script 48
#
# Usage: bash src/run_downstream_chain.sh

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
echo "Downstream Dependency Chain"
echo "============================================"

# ── Step 0: Collect all existing chr4-22 COLOC job IDs ──
EXISTING_COLOC_JOBS=$(squeue --me -n susie_coloc --format="%i" 2>/dev/null | grep -o "^[0-9]*" | sort -un | tr '\n' ',' | sed 's/,$//')
echo "Existing chr4-22 COLOC jobs: ${EXISTING_COLOC_JOBS}"

# ── Step 1: Wait for eQTL SuSiE to finish, then submit chr1-3 COLOC ──
echo ""
echo "--- Waiting for eQTL SuSiE chr1-3 to finish ---"
ALL_CHR13_JOBS=""

for chr in 1 2 3; do
  TOTAL_GENES=0
  case $chr in
    1) TOTAL_GENES=2048 ;;
    2) TOTAL_GENES=1284 ;;
    3) TOTAL_GENES=1090 ;;
  esac

  while true; do
    count=$(find results/eqtl_susie/chr${chr}/ -name "*.rds" 2>/dev/null | wc -l)
    workers=$(squeue --me -n eqtl_susie --format="%i" 2>/dev/null | grep "_${chr}$" | wc -l)

    if [ $workers -eq 0 ] || [ $count -ge $TOTAL_GENES ]; then
      echo "  chr${chr}: READY (${count}/${TOTAL_GENES} RDS, ${workers} workers remaining)"
      # Memory: chr1,2 need 128G (large LD blocks); chr3 is safe at 96G
      if [ $chr -le 2 ]; then
        CHR_MEM="128G"
      else
        CHR_MEM="96G"
      fi
      echo "  Submitting chr${chr} COLOC at ${CHR_MEM}"
      CHR_JOBS=""
      for gwas in "${GWAS_LIST[@]}"; do
        JID=$(GWAS_NAME="${gwas}" sbatch --parsable --mem=${CHR_MEM} --array=${chr} src/06_susie_coloc.sh)
        CHR_JOBS="${CHR_JOBS}${CHR_JOBS:+,}${JID}"
      done
      ALL_CHR13_JOBS="${ALL_CHR13_JOBS}${ALL_CHR13_JOBS:+,}${CHR_JOBS}"
      echo "  Submitted chr${chr} COLOC: ${#GWAS_LIST[@]} jobs at ${CHR_MEM}"
      break
    fi

    echo "  chr${chr}: ${count}/${TOTAL_GENES} (${workers} workers) — waiting 2 min..."
    sleep 120
  done
done

echo ""
echo "All chr1-3 COLOC submitted."

# ── Step 2: Build dependency string for ALL COLOC jobs ──
ALL_COLOC="${EXISTING_COLOC_JOBS}${EXISTING_COLOC_JOBS:+,}${ALL_CHR13_JOBS}"
# Convert to afterok format
DEP_STR=$(echo "$ALL_COLOC" | tr ',' '\n' | sort -un | tr '\n' ':' | sed 's/:$//')
echo "Dependency chain on ${DEP_STR}"

# ── Step 3: Script 07 — combine all COLOC results ──
echo ""
echo "--- Submitting Script 07 (combine COLOC) ---"
JOB_07=$(sbatch --parsable \
  --dependency=afterok:${DEP_STR} \
  --partition=cpu --mem=32G --time=4:00:00 \
  --job-name=combine_coloc \
  --output=logs/combine_coloc_%j.out \
  --error=logs/combine_coloc_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate finemapping && Rscript src/07_combine_susie_coloc.R")
echo "  Job 07 (combine COLOC): ${JOB_07} (afterok: all COLOC)"

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
echo "  Job 27a (atlas): ${JOB_27a} (afterok: ${JOB_07})"

# ── Step 5: Script 75 — integrate causal overhaul ──
echo ""
echo "--- Submitting Script 75 (causal overhaul integration) ---"
JOB_75=$(sbatch --parsable \
  --dependency=afterok:${JOB_27a} \
  --partition=cpu --mem=64G --time=8:00:00 \
  --job-name=causal_75 \
  --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/causal_75_%j.out \
  --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/causal_75_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate rnaseq && Rscript RNA-seq/75_integrate_causal_overhaul.R")
echo "  Job 75 (causal): ${JOB_75} (afterok: ${JOB_27a})"

# ── Step 6: Script 48 — cross-ancestry replication ──
echo ""
echo "--- Submitting Script 48 (cross-ancestry replication) ---"
JOB_48=$(sbatch --parsable \
  --dependency=afterok:${JOB_75} \
  --partition=cpu --mem=32G --time=4:00:00 \
  --job-name=cross_ancestry_48 \
  --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cross_ancestry_48_%j.out \
  --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cross_ancestry_48_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate rnaseq && Rscript RNA-seq/48_cross_ancestry_replication.R")
echo "  Job 48 (cross-ancestry): ${JOB_48} (afterok: ${JOB_75})"

# ── Step 7: Script 217 — stratified causal atlas integration ──
echo ""
echo "--- Submitting Script 217 (stratified causal atlas) ---"
JOB_217=$(sbatch --parsable \
  --dependency=afterok:${JOB_75} \
  --partition=cpu --mem=32G --time=4:00:00 \
  --job-name=stratified_217 \
  --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/stratified_217_%j.out \
  --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/stratified_217_%j.err \
  --wrap="cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design && eval \"\$(micromamba shell hook -s bash)\" && micromamba activate rnaseq && Rscript RNA-seq/217_stratified_causal_atlas.R")
echo "  Job 217 (stratified): ${JOB_217} (afterok: ${JOB_75})"

echo ""
echo "============================================"
echo "FULL DEPENDENCY CHAIN:"
echo "  All COLOC (chr4-22 + chr1-3)"
echo "    └→ Script 07 (combine)     : ${JOB_07}"
echo "        └→ Script 27a (atlas)  : ${JOB_27a}"
echo "            └→ Script 75       : ${JOB_75}"
echo "                ├→ Script 48   : ${JOB_48}"
echo "                └→ Script 217  : ${JOB_217}"
echo ""
echo "Monitor: squeue --me"
echo "============================================"
