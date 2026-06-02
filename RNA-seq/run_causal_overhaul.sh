#!/bin/bash
# =============================================================================
# MR / TWAS / COLOC Pipeline Overhaul — Master Orchestrator
#
# Submits all phases with proper SLURM dependency chains.
# Designed to be idempotent: already-completed steps are skipped.
#
# Usage:
#   bash RNA-seq/run_causal_overhaul.sh
#
# Monitor:
#   squeue -u $USER --format="%.18i %.25j %.8T %.10M"
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

mkdir -p RNA-seq/logs

echo "============================================================"
echo "  MR/TWAS/COLOC Pipeline Overhaul"
echo "  Start: $(date)"
echo "============================================================"

# Track all job IDs for final summary
declare -A JOBS

# ============================================================================
# Phase 0D: GWAS Liftover (if not done)
# ============================================================================
HG19_DIR="${BASE}/GWAS/MR_Data/hg19"
if [[ ! -f "${HG19_DIR}/Ghodsian_2021_NAFLD_harmonised_hg19.tsv.gz" ]]; then
  echo ""
  echo "--- Phase 0D: GWAS Liftover ---"
  JOB=$(sbatch --parsable RNA-seq/sbatch_liftover_gwas.sbatch)
  JOBS[liftover]=${JOB}
  echo "  Job: ${JOB}"
  LIFTOVER_DEP="--dependency=afterok:${JOB}"
else
  echo "Phase 0D: GWAS liftover already complete"
  LIFTOVER_DEP=""
fi

# ============================================================================
# Phase 1A: Enhanced cis-MR (Ghodsian + UKBB enzymes)
# ============================================================================
echo ""
echo "--- Phase 1A: Enhanced MR ---"

for GWAS in ghodsian ukbb_alt ukbb_ast ukbb_ggt pdff; do
  JOB=$(sbatch --parsable ${LIFTOVER_DEP} \
    --export=ALL,GWAS_NAME="${GWAS}" \
    --job-name="mr_${GWAS}" \
    --output="RNA-seq/logs/enhanced_mr_${GWAS}_%j.out" \
    --error="RNA-seq/logs/enhanced_mr_${GWAS}_%j.err" \
    RNA-seq/run_enhanced_mr.sbatch)
  JOBS[mr_${GWAS}]=${JOB}
  echo "  ${GWAS}: Job ${JOB}"
done

# ============================================================================
# Phase 1B: Bidirectional MR (depends on liftover)
# ============================================================================
echo ""
echo "--- Phase 1B: Bidirectional MR ---"
JOB=$(sbatch --parsable ${LIFTOVER_DEP} RNA-seq/run_bidirectional_mr.sbatch)
JOBS[bidir_mr]=${JOB}
echo "  Job: ${JOB}"

# ============================================================================
# Phase 2A-Step1: OTTERS Format Broadaway (no dependencies)
# ============================================================================
echo ""
echo "--- Phase 2A: OTTERS Format ---"
OTTERS_FMT_DIR="${BASE}/data/broadaway_eqtl/otters_format"
if [[ ! -f "${OTTERS_FMT_DIR}/chr1_broadaway.txt" ]]; then
  JOB=$(sbatch --parsable RNA-seq/run_otters_format.sbatch)
  JOBS[otters_fmt]=${JOB}
  echo "  Job: ${JOB}"
  OTTERS_FMT_DEP="--dependency=afterok:${JOB}"
else
  echo "  OTTERS format already complete"
  OTTERS_FMT_DEP=""
fi

# ============================================================================
# Phase 3A: SuSiE-COLOC Rerun (no dependency — uses existing Broadaway data)
# ============================================================================
echo ""
echo "--- Phase 3A: SuSiE-COLOC Rerun ---"
GWAS_LIST=(UKBB_ALT UKBB_AST UKBB_GGT FINNGEN_NAFLD FINNGEN_NASH FINNGEN_HCC)
SUSIE_JOBS=""

for GWAS_NAME in "${GWAS_LIST[@]}"; do
  JOB=$(sbatch --parsable \
    --job-name="susie_${GWAS_NAME}" \
    --output="RNA-seq/logs/susie_${GWAS_NAME}_%j.out" \
    --error="RNA-seq/logs/susie_${GWAS_NAME}_%j.err" \
    --export=ALL,GWAS_NAME="${GWAS_NAME}" \
    RNA-seq/sbatch_susie_coloc.sh)
  JOBS[susie_${GWAS_NAME}]=${JOB}
  SUSIE_JOBS="${SUSIE_JOBS}:${JOB}"
  echo "  ${GWAS_NAME}: Job ${JOB}"
done

# ============================================================================
# Phase 2A-Step2: OTTERS Weight Training (depends on format)
# ============================================================================
echo ""
echo "--- Phase 2A: OTTERS Weight Training ---"
OTTERS_WEIGHT_DIR="${BASE}/data/broadaway_eqtl/otters_weights"
if [[ ! -d "${OTTERS_WEIGHT_DIR}/chr22" ]]; then
  JOB=$(sbatch --parsable ${OTTERS_FMT_DEP} RNA-seq/run_otters_train.sbatch)
  JOBS[otters_train]=${JOB}
  echo "  Job: ${JOB} (array 1-22)"
  OTTERS_TRAIN_DEP="--dependency=afterok:${JOB}"
else
  echo "  OTTERS weights already complete"
  OTTERS_TRAIN_DEP=""
fi

# ============================================================================
# Phase 2A-Step3: OTTERS TWAS (depends on training + liftover)
# ============================================================================
echo ""
echo "--- Phase 2A: OTTERS TWAS ---"
COMBINED_DEP=""
if [[ -n "${OTTERS_TRAIN_DEP}" ]] || [[ -n "${LIFTOVER_DEP}" ]]; then
  DEP_LIST=""
  [[ -n "${OTTERS_TRAIN_DEP}" ]] && DEP_LIST="${OTTERS_TRAIN_DEP#--dependency=afterok:}"
  [[ -n "${LIFTOVER_DEP}" ]] && DEP_LIST="${DEP_LIST:+${DEP_LIST}:}${LIFTOVER_DEP#--dependency=afterok:}"
  [[ -n "${DEP_LIST}" ]] && COMBINED_DEP="--dependency=afterok:${DEP_LIST}"
fi

for GWAS in ghodsian ukbb_alt ukbb_ast ukbb_ggt; do
  JOB=$(sbatch --parsable ${COMBINED_DEP} \
    --export=ALL,GWAS_NAME="${GWAS}" \
    --job-name="otwas_${GWAS}" \
    --output="RNA-seq/logs/otters_twas_${GWAS}_%j.out" \
    --error="RNA-seq/logs/otters_twas_${GWAS}_%j.err" \
    RNA-seq/run_otters_twas.sbatch)
  JOBS[otters_twas_${GWAS}]=${JOB}
  echo "  ${GWAS}: Job ${JOB}"
done

# ============================================================================
# Phase 2B: cTWAS with GTEx Liver weights (depends on liftover)
# ============================================================================
echo ""
echo "--- Phase 2B: cTWAS ---"
for GWAS in ghodsian ukbb_alt; do
  JOB=$(sbatch --parsable ${LIFTOVER_DEP} \
    --export=ALL,GWAS_NAME="${GWAS}" \
    --job-name="ctwas_${GWAS}" \
    --output="RNA-seq/logs/ctwas_${GWAS}_%j.out" \
    --error="RNA-seq/logs/ctwas_${GWAS}_%j.err" \
    RNA-seq/run_ctwas.sbatch)
  JOBS[ctwas_${GWAS}]=${JOB}
  echo "  ${GWAS}: Job ${JOB}"
done

# ============================================================================
# Phase 3B: Multi-Trait COLOC (depends on liftover for hg19 GWAS)
# ============================================================================
echo ""
echo "--- Phase 3B: Multi-Trait COLOC ---"
JOB=$(sbatch --parsable ${LIFTOVER_DEP} RNA-seq/run_hyprcoloc.sbatch)
JOBS[hyprcoloc]=${JOB}
echo "  Job: ${JOB}"

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "============================================================"
echo "  All jobs submitted!"
echo "============================================================"
echo ""
echo "Job summary:"
for key in "${!JOBS[@]}"; do
  printf "  %-20s : %s\n" "${key}" "${JOBS[$key]}"
done

echo ""
echo "Monitor: squeue -u \$USER --format='%.18i %.25j %.8T %.10M'"
echo ""
echo "After all complete, run Phase 5 (atlas integration):"
echo "  Rscript RNA-seq/75_integrate_causal_overhaul.R"
echo ""
