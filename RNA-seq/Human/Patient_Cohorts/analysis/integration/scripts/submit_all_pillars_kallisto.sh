#!/usr/bin/env bash
# submit_all_pillars_kallisto.sh
# ---------------------------------------------------------------------------
# Master submission for the 4-pillar DEG robustness framework on kallisto
# dream results. Submits all SLURM jobs with proper dependencies.
#
# Uses standalone sbatch scripts (NOT --wrap) to ensure #!/bin/bash shebang
# is honored. The --wrap approach defaults to /bin/sh which breaks
# `set -eo pipefail`.
#
# Prerequisites:
#   - merged_dge.rds has been rebuilt from kallisto counts (build_kallisto_dge.R)
#   - dream_results.csv is the canonical kallisto dream output (39,276 genes)
#   - STAR iteration files have been moved to *_star_backup/ dirs
# ---------------------------------------------------------------------------
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGDIR=${SCRIPTS}/logs
mkdir -p "${LOGDIR}"

cd "${SCRIPTS}"

echo "=================================================="
echo " 4-Pillar DEG Robustness — Kallisto Rebuild"
echo " $(date)"
echo "=================================================="

# ============================================================================
# Pillar A: CPSS (100 pairs = 200 dream fits)
# ============================================================================
echo ""
echo "--- Pillar A: CPSS (100 pairs) ---"
CPSS_JID=$(sbatch --parsable run_cpss.sbatch)
echo "  CPSS array: ${CPSS_JID}"

# ============================================================================
# Pillar A: Bootstrap (1000 iters)
# ============================================================================
echo ""
echo "--- Pillar A: Bootstrap (1000 iters) ---"
BOOT_JID=$(sbatch --parsable run_bootstrap.sbatch)
echo "  Bootstrap array: ${BOOT_JID}"

# ============================================================================
# Pillar A: K-fold (60 iters: 10 for K=10, 50 for K=50)
# ============================================================================
echo ""
echo "--- Pillar A: K-fold (60 iters) ---"
KFOLD_JID=$(sbatch --parsable run_kfold.sbatch)
echo "  K-fold array: ${KFOLD_JID}"

# ============================================================================
# Pillar A: LTO (10 pairs from 5 cohorts)
# ============================================================================
echo ""
echo "--- Pillar A: LTO (10 pairs) ---"
LTO_JID=$(sbatch --parsable run_lto.sbatch)
echo "  LTO array: ${LTO_JID}"

# ============================================================================
# LOO-CV (5 folds — needed for Pillar B)
# ============================================================================
echo ""
echo "--- LOO-CV (5 folds, prerequisite for Pillar B) ---"
COHORTS=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)
LOO_JIDS=()
for cohort in "${COHORTS[@]}"; do
    JID=$(sbatch --parsable --export=ALL,HELD_OUT="${cohort}" run_loco.sbatch)
    echo "  LOO ${cohort}: ${JID}"
    LOO_JIDS+=("${JID}")
done

# ============================================================================
# Pillar B: Transferability (depends on LOO-CV)
# ============================================================================
echo ""
echo "--- Pillar B: LOCO transferability (depends on LOO-CV) ---"
LOO_DEP=$(IFS=:; echo "${LOO_JIDS[*]}")
PILLAR_B_JID=$(sbatch --parsable --dependency=afterok:${LOO_DEP} run_pillarB.sbatch)
echo "  Pillar B: ${PILLAR_B_JID} (depends on LOO: ${LOO_DEP})"

# ============================================================================
# Pillar C: Permutation null (1000 iters)
# ============================================================================
echo ""
echo "--- Pillar C: Permutation null (1000 iters) ---"
PERM_JID=$(sbatch --parsable run_permutation.sbatch)
echo "  Permutation array: ${PERM_JID}"

# ============================================================================
# Pillar D: PVCA
# ============================================================================
echo ""
echo "--- Pillar D: PVCA ---"
PVCA_JID=$(sbatch --parsable run_pvca.sbatch)
echo "  PVCA: ${PVCA_JID}"

# ============================================================================
# Pillar D: SVA sensitivity
# ============================================================================
echo ""
echo "--- Pillar D: SVA ---"
SVA_JID=$(sbatch --parsable run_sva.sbatch)
echo "  SVA: ${SVA_JID}"

# ============================================================================
# Pillar D: SVA v2 (residuals-based)
# ============================================================================
echo ""
echo "--- Pillar D: SVA v2 ---"
SVA2_JID=$(sbatch --parsable run_sva_v2.sbatch)
echo "  SVA v2: ${SVA2_JID}"

# ============================================================================
# Aggregation jobs (depend on iteration arrays)
# ============================================================================
echo ""
echo "=== Aggregation (chained dependencies) ==="

# Pillar A aggregation (after CPSS + bootstrap + kfold)
echo ""
echo "--- Pillar A: Aggregation (after CPSS/bootstrap/kfold) ---"
AGG_A_JID=$(sbatch --parsable --dependency=afterok:${CPSS_JID}:${BOOT_JID}:${KFOLD_JID} run_pillarA_agg.sbatch)
echo "  Pillar A agg: ${AGG_A_JID} (depends on CPSS:${CPSS_JID} + boot:${BOOT_JID} + kfold:${KFOLD_JID})"

# Pillar A LTO aggregation (after LTO)
echo ""
echo "--- Pillar A: LTO aggregation (after LTO) ---"
AGG_LTO_JID=$(sbatch --parsable --dependency=afterok:${LTO_JID} run_pillarA_lto_agg.sbatch)
echo "  LTO agg: ${AGG_LTO_JID} (depends on LTO:${LTO_JID})"

# Pillar C aggregation (after permutation)
echo ""
echo "--- Pillar C: Aggregation (after permutation) ---"
AGG_C_JID=$(sbatch --parsable --dependency=afterok:${PERM_JID} run_pillarC_agg.sbatch)
echo "  Pillar C agg: ${AGG_C_JID} (depends on perm:${PERM_JID})"

echo ""
echo "=================================================="
echo " All jobs submitted."
echo ""
echo " Job IDs:"
echo "   CPSS array:        ${CPSS_JID}"
echo "   Bootstrap array:   ${BOOT_JID}"
echo "   K-fold array:      ${KFOLD_JID}"
echo "   LTO array:         ${LTO_JID}"
echo "   LOO-CV:            ${LOO_JIDS[*]}"
echo "   Pillar B:          ${PILLAR_B_JID}"
echo "   Permutation array: ${PERM_JID}"
echo "   PVCA:              ${PVCA_JID}"
echo "   SVA:               ${SVA_JID}"
echo "   SVA v2:            ${SVA2_JID}"
echo "   Pillar A agg:      ${AGG_A_JID}"
echo "   LTO agg:           ${AGG_LTO_JID}"
echo "   Pillar C agg:      ${AGG_C_JID}"
echo ""
echo " Monitor: squeue -u \$(whoami) --name=CPSS,bootstrap,k-fold,LTO,LOCO,permutation-null,PVCA,SVA"
echo "=================================================="
