#!/bin/bash
#SBATCH --job-name=ComBatseq
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/34c2_combat_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/34c2_combat_%j.err

set -eo pipefail
source ~/.bashrc
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-8}
Rscript RNA-seq/34c2_combat_seq_sensitivity_lvqw.R
