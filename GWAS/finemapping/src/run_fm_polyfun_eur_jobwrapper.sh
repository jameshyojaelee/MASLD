#!/bin/bash -l
# Owner per-locus job wrapper for Phase 9e PolyFun standalone FM.
# need micromamba initialized in their shell).
#
# Identical job logic to src/03_run_fm_per_locus.sh, but uses absolute paths
# micromamba install via absolute paths so any user with read access to
# /gpfs/commons/home/jameslee can run it without environment setup.
#
# Called by run_fm_polyfun_eur.sh as an sbatch script.

set -o pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping

# Initialize micromamba via absolute path (works for any user with
# read access; bypasses need for micromamba in their PATH)
export MAMBA_EXE="/gpfs/commons/home/jameslee/.local/bin/micromamba"
export MAMBA_ROOT_PREFIX="/gpfs/commons/home/jameslee/micromamba"
eval "$(${MAMBA_EXE} shell hook -s bash)"
micromamba activate finemapping

if [ -z "${CONDA_PREFIX:-}" ]; then
  echo "ERROR: micromamba activation failed. Check that ${MAMBA_EXE} is readable."
  exit 1
fi

set -u

sumstats_name=$1
ld_pop=$2
lead_snps_file=$3
N_tot=$4
N_cases=$5
window_mb=$6
ancestry=$7

LOCUS=$(awk -v line=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==line {print $3}' "${lead_snps_file}")

if [ -z "$LOCUS" ]; then
  echo "ERROR: No locus found for array task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "[polyfun-fm] task=${SLURM_ARRAY_TASK_ID}  study=${sumstats_name}  ld=${ld_pop}  locus=${LOCUS}  ancestry=${ancestry}"
echo "[polyfun-fm] env=$(basename ${CONDA_PREFIX})  user=$(whoami)  start=$(date)"

Rscript src/03_run_fm_per_locus.R \
    "${sumstats_name}" \
    "${ld_pop}" \
    "${LOCUS}" \
    "${N_tot}" \
    "${N_cases}" \
    "${window_mb}" \
    "${ancestry}"

rc=$?
echo "[polyfun-fm] task=${SLURM_ARRAY_TASK_ID}  locus=${LOCUS}  done=$(date)  rc=${rc}"
exit ${rc}
