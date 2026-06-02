#!/usr/bin/env bash
# =====================================================================
# run_mash_loo_sweep.sh  (packed bigmem version, 3 stages)
#
# Submits 3 packed SLURM jobs that share one whole bigmem node each
# (72 CPU + 512G RAM, clears MinTRES=mem=500G). The actual stage bodies
# live in sibling scripts so they run as proper bash, not /bin/sh.
#
#   Stage 1: _stage1_megas.sh   (no deps)
#   Stage 2: _stage2_loos.sh    (no deps; parallel with Stage 1)
#   Stage 3: _stage3_sweeps.sh  (afterok:Stage 1 + Stage 2)
#
# Stage 2 does NOT depend on Stage 1: dream_loo_cv_contrasts.R rebuilds
# DGE from merged_counts_raw.rds and does not consume Stage 1 output.
# Stage 3 consumes both.
# =====================================================================

set -euo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
INT=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration
SCRIPTS=${INT}/scripts
LOGS=${INT}/logs/mash_loo_sweep
mkdir -p "${LOGS}"

COMMON_NODE="--partition=bigmem --cpus-per-task=72 --mem=512G --time=48:00:00 --nodes=1 --exclusive"
COMMON_SMALL="--partition=bigmem --cpus-per-task=16 --mem=512G --time=48:00:00 --nodes=1"

echo "============================================================"
echo "MASH LOO + LFC sweep orchestrator  (packed bigmem)"
echo "Submitted: $(date)"
echo "Logs dir: ${LOGS}"
echo "============================================================"

# Stage 1
stage1_id=$(sbatch --parsable \
    --job-name=mash_stage1_megas \
    --output="${LOGS}/stage1_%j.out" \
    --error="${LOGS}/stage1_%j.err" \
    ${COMMON_NODE} \
    "${SCRIPTS}/_stage1_megas.sh")
echo "  -> ${stage1_id} (stage1_megas, 72c/512G/excl, no dep)"

# Stage 2
stage2_id=$(sbatch --parsable \
    --job-name=mash_stage2_loos \
    --output="${LOGS}/stage2_%j.out" \
    --error="${LOGS}/stage2_%j.err" \
    ${COMMON_NODE} \
    "${SCRIPTS}/_stage2_loos.sh")
echo "  -> ${stage2_id} (stage2_loos,  72c/512G/excl, no dep)"

# Stage 3
stage3_id=$(sbatch --parsable \
    --job-name=mash_stage3_sweeps \
    --output="${LOGS}/stage3_%j.out" \
    --error="${LOGS}/stage3_%j.err" \
    --dependency=afterok:${stage1_id}:${stage2_id} \
    ${COMMON_SMALL} \
    "${SCRIPTS}/_stage3_sweeps.sh")
echo "  -> ${stage3_id} (stage3_sweeps, 16c/512G, dep ${stage1_id}+${stage2_id})"

echo ""
echo "============================================================"
echo "Job IDs:"
printf "  stage1_megas   = %s\n" "${stage1_id}"
printf "  stage2_loos    = %s\n" "${stage2_id}"
printf "  stage3_sweeps  = %s\n" "${stage3_id}"
echo "============================================================"
echo "Monitor: squeue -u \$USER -o '%.10i %.22j %.2t %.10M %.6D %R'"
echo "Logs:    ${LOGS}"
