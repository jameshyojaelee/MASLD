#!/bin/bash
#SBATCH --job-name=v6_P1
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --output=logs/P1_M2_random_%j.out
#SBATCH --error=logs/P1_M2_random_%j.err
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
cd "${SCRIPT_DIR}"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "=== P1 dream M2 random-slope refit start: $(date) ==="
Rscript 10_M2_refit_random.R
echo "=== done: $(date) ==="
