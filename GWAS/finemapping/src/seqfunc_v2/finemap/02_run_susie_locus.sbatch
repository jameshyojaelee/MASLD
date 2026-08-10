#!/bin/bash -l
#SBATCH --job-name=SuSiE
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
# Log paths are RELATIVE and resolved against the submission working
# directory, which submit_uniform35_v2.sh sets to $UNIFORM35_RUN_ROOT.
# They were absolute paths into runs/uniform35_v2_2026-07-13/logs/, which
# SLURM evaluates before the script body exports the env var -- so a v3 run
# silently wrote its logs into the v2 run directory.  The submit script also
# passes --output/--error explicitly, which overrides these.
#SBATCH --output=logs/susie_%A_%a.out
#SBATCH --error=logs/susie_%A_%a.err

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
FM_ROOT="${PROJECT_ROOT}/GWAS/finemapping"
RUN_ROOT="${UNIFORM35_RUN_ROOT:-${FM_ROOT}/runs/uniform35_v2_2026-07-13}"
TASK_FILE="${1:?Usage: 02_run_susie_locus.sbatch <task_file>}"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export UNIFORM35_RUN_ROOT="${RUN_ROOT}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

LOCUS_INDEX=$(awk -F '\t' -v row="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == row {print $2}' "${TASK_FILE}")
if [[ -z "${LOCUS_INDEX}" ]]; then
  echo "No locus_index for task ${SLURM_ARRAY_TASK_ID} in ${TASK_FILE}" >&2
  exit 2
fi

eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping
exec Rscript "${FM_ROOT}/src/seqfunc_v2/finemap/02_run_susie_locus.R" "${LOCUS_INDEX}"
