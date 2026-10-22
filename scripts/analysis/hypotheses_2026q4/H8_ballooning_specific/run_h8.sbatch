#!/bin/bash
#SBATCH --job-name=h8_ballooning
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H8-20261009T015529Z/slurm-%x-%j.out
# Usage: sbatch run_h8.sbatch <label> <n_perm>
# Writes Analysis/Hypotheses/H8-20261009T015529Z/<label>-<UTC timestamp>/ (new directory).
set -euo pipefail
label=${1:?label}
n_perm=${2:?n_perm}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
rscript=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
hyp=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H8-20261009T015529Z
script=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H8_ballooning_specific/h8_ballooning_meta.R
out="${hyp}/${label}-$(date -u +%Y%m%dT%H%M%SZ)"
echo "host $(hostname) job ${SLURM_JOB_ID} out ${out}"
"${rscript}" "${script}" "${out}" "${n_perm}" "${SLURM_CPUS_PER_TASK}"
cp "${script}" "${out}/h8_ballooning_meta.as_run.R"
sha256sum "${script}" > "${out}/script.sha256"
