#!/bin/bash
#SBATCH --job-name=dream-bariatric
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto/logs/06a_bariatric_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto/logs/06a_bariatric_%j.err

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

echo "[$(date)] Starting bariatric sensitivity arm (06a)"
echo "Host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Mem: ${SLURM_MEM_PER_NODE}MB"

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto/06a_bariatric_sensitivity.R

echo "[$(date)] Done"
