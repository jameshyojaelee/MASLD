#!/bin/bash
#SBATCH --job-name=snapatac2
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=logs/zonation_F064_rerun_%j.out
#SBATCH --error=logs/zonation_F064_rerun_%j.err
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome
mkdir -p logs
echo "==== F064 re-run: zonation-resolved hepatocyte DA (script 15) only ===="
python scripts/15_zonation_da_hepatocytes.py
echo "ALL DONE"
