#!/bin/bash
#SBATCH --job-name=SuSiE
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=logs/one_%j.out
#SBATCH --error=logs/one_%j.err

# Run a single locus by index.  A real script rather than `sbatch --wrap`,
# which executes under /bin/sh here and cannot parse the bash-only micromamba
# shell hook.
#
#   sbatch --export=ALL 02_run_one.sbatch <locus_index>
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
RUN_ROOT="${UNIFORM35_RUN_ROOT:?UNIFORM35_RUN_ROOT must be set}"
export MASLD_PROJECT_ROOT UNIFORM35_RUN_ROOT

export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

LOCUS_INDEX="${1:?usage: 02_run_one.sbatch <locus_index>}"
echo "[$(date -u +%FT%TZ)] single locus ${LOCUS_INDEX} in ${RUN_ROOT}"
Rscript "${PROJECT_ROOT}/GWAS/finemapping/src/seqfunc_v2/finemap/02_run_susie_locus.R" "${LOCUS_INDEX}"
