#!/bin/bash
#SBATCH --job-name=s_v3_M08
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=logs/M08_%j.out
#SBATCH --error=logs/M08_%j.err
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
cd "${SCRIPT_DIR}"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "=== Module 08 (stratified-ashr) start: $(date) ==="
Rscript 08_stratified_ashr.R
echo "=== Module 08 end: $(date) ==="
