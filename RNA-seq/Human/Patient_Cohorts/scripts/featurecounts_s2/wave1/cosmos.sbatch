#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --job-name=COSMOS
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/COSMOS_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/COSMOS_%j.err
# Wave1 [C] TERMINAL: 52 COSMOS mechanistic paths. Reads -s2 dream. Env = rnaseq (script self-contained).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 52 COSMOS mechanistic paths | $(date) ==="
Rscript RNA-seq/52_cosmos_mechanistic_paths.R
if [ $? -ne 0 ]; then echo "FAILED: 52_cosmos_mechanistic_paths"; exit 1; fi
echo "=== done $(date) ==="
