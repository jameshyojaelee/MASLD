#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/loo_s2_2026-05-29/dream_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/loo_s2_2026-05-29/dream_%A_%a.err
# =====================================================================
# run_loo_fold_array.sh  (2026-05-29, deadlock-safe STAR -s 2 rebuild)
#
# Resubmission of the mash/masl LOO-CV folds AFTER run_mash_loo_sweep.sh
# deadlocked (it backgrounded ~7 concurrent dream() calls in ONE bigmem
# node -> BiocParallel/OpenBLAS thread exhaustion). THE FIX: exactly ONE
# dream() per SLURM allocation. Here that is one dream() per ARRAY TASK
# (each task = its own node allocation). No '&' / backgrounding anywhere.
#
# Each task reads one (CONTRAST, MASH_DEF, HELD_OUT) row from folds.tsv
# (1-based data rows -> array index 0..N-1) and runs the canonical
# per-fold producer dream_loo_cv_contrasts.R, which reads the canonical
# -s 2 merged_counts_raw.rds / meta_matched.rds and writes
# results/integration/loo_cv/<subdir>/dream_loo_<cohort>.csv.
#
# BLAS thread cap (OMP/OPENBLAS=8 + SLURM_CPUS_PER_TASK=8) prevents the
# RLIMIT_NPROC blowup that compounds the MulticoreParam stall.
# =====================================================================
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
INT=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration
SCRIPTS=${INT}/scripts
MANIFEST=${SCRIPTS}/loo_s2_2026-05-29/folds.tsv

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

# --- select this task's fold from the manifest (skip header) ---
row=$(( SLURM_ARRAY_TASK_ID + 2 ))   # +1 for header, +1 for 1-based sed
line=$(sed -n "${row}p" "${MANIFEST}")
if [[ -z "${line}" ]]; then
    echo "ERROR: no manifest row for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} (line ${row})" >&2
    exit 1
fi
CONTRAST=$(echo "${line}" | cut -f1)
MASH_DEF_RAW=$(echo "${line}" | cut -f2)
HELD_OUT=$(echo "${line}" | cut -f3)

# masl_vs_healthy ignores MASH_DEF but the producer validates it against
# {borderline_grouped, strict}; map the manifest sentinel 'none' to the
# producer default. (For masl_vs_healthy the value is never read.)
if [[ "${MASH_DEF_RAW}" == "none" ]]; then
    MASH_DEF="borderline_grouped"
else
    MASH_DEF="${MASH_DEF_RAW}"
fi

# --- cap BLAS / OpenMP threads to the allocation (8) ---
export SLURM_CPUS_PER_TASK=8
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export R_PARALLEL_SEED=42

export CONTRAST MASH_DEF HELD_OUT

echo "=== LOO fold (STAR -s 2) ==="
echo "host=$(hostname)  date=$(date)"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}  SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}  JOB=${SLURM_JOB_ID}"
echo "CONTRAST=${CONTRAST}  MASH_DEF=${MASH_DEF} (raw=${MASH_DEF_RAW})  HELD_OUT=${HELD_OUT}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}  OMP_NUM_THREADS=${OMP_NUM_THREADS}  OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
echo ""

cd "${PROJECT}"

# EXACTLY ONE dream() per job. No backgrounding.
Rscript "${SCRIPTS}/dream_loo_cv_contrasts.R"
rc=$?

echo ""
echo "=== fold finished rc=${rc}  date=$(date) ==="
exit ${rc}
