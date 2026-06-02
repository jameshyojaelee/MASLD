#!/usr/bin/env bash
#SBATCH --job-name=bayev46d
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/46d_%j.out
#SBATCH --error=RNA-seq/logs/46d_%j.err

set -eo pipefail   # NB: no -u; conda/micromamba activate scripts reference
                   # unbound variables in their standard setup
source ~/.bashrc
micromamba activate rnaseq

cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p RNA-seq/logs

# Full permutation count for production
export N_PERM="${N_PERM:-1000}"

echo "=== 46d Convergence Evidence — SLURM job $SLURM_JOB_ID ==="
echo "Date: $(date)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "N_PERM: $N_PERM"

time Rscript RNA-seq/46d_convergence_evidence.R

echo "=== Done ==="
