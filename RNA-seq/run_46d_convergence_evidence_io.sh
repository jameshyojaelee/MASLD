#!/usr/bin/env bash
#SBATCH --job-name=46d-conv-io
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=RNA-seq/logs/46d_io_%j.out
#SBATCH --error=RNA-seq/logs/46d_io_%j.err

set -eo pipefail

source ~/.bashrc
set +u
micromamba activate rnaseq
set -u

cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p RNA-seq/logs

export N_PERM="${N_PERM:-1000}"

echo "=== 46d Convergence Evidence (io) — SLURM job $SLURM_JOB_ID ==="
echo "Date: $(date)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "N_PERM: $N_PERM"

time Rscript RNA-seq/46d_convergence_evidence.R

echo "=== Done at $(date) ==="
