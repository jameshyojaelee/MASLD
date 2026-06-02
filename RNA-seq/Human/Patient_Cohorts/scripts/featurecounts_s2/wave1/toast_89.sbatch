#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --job-name=TOAST
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/TOAST_89_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/TOAST_89_%j.err
# Wave1 [F]: 89 TOAST cell-type-resolved DE. Reads -s2 merged_dge.rds. Internal per-CT loop.
# NOT an atlas feeder (27a/75/217 do not read toast_celltype_de_all.csv).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 89 TOAST celltype DE | $(date) ==="
Rscript RNA-seq/89_toast_celltype_de.R
if [ $? -ne 0 ]; then echo "FAILED: 89_toast_celltype_de"; exit 1; fi
echo "=== done $(date) ==="
