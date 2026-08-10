#!/bin/bash
#SBATCH --job-name=SuSiE
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pilot_%A_%a.out
#SBATCH --error=logs/pilot_%A_%a.err

# Runs one pilot locus per array task.  Memory is the xlarge tier for every task
# because the pilot deliberately spans the whole size range and it is only 14
# jobs; the full array uses the per-tier request instead.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
# MUST match 02_run_susie_locus.sbatch (finemapping, data.table 1.14.10).
# rnaseq carries data.table 1.18.2.1, whose zero-length-column recycling
# differs: data.table(x=character(0), y="a") yields 1 NA-filled row under
# 1.14.10 and 0 rows under 1.18.2. Running the same script under a different
# env silently changes results.
micromamba activate finemapping
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SRC="${PROJECT_ROOT}/GWAS/finemapping/src/seqfunc_v2/finemap"
RUN_ROOT="${UNIFORM35_RUN_ROOT:?UNIFORM35_RUN_ROOT must be set}"
export MASLD_PROJECT_ROOT UNIFORM35_RUN_ROOT

# Threads for the BLAS-3 Cholesky ladder; SuSiE itself stays effectively serial.
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

PILOT="${RUN_ROOT}/config/pilot_loci.tsv"
LOCUS_INDEX=$(awk -F '\t' -v row="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == row {print $1}' "${PILOT}")
[[ -n "${LOCUS_INDEX}" ]] || { echo "no locus for array task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }

echo "[$(date -u +%FT%TZ)] pilot task ${SLURM_ARRAY_TASK_ID} -> locus_index ${LOCUS_INDEX}"
Rscript "${SRC}/02_run_susie_locus.R" "${LOCUS_INDEX}"
