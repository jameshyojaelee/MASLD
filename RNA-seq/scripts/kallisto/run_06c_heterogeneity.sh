#!/bin/bash
#SBATCH --job-name=metafor-isq
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto/logs/06c_heterogeneity_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto/logs/06c_heterogeneity_%j.err

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

echo "[$(date)] Starting I2 heterogeneity investigation (06c)"
echo "Host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Mem: ${SLURM_MEM_PER_NODE}MB"

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto/06c_heterogeneity_investigation.R

echo "[$(date)] Done"
