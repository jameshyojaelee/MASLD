#!/bin/bash
#SBATCH --job-name=wave0_dream_kallisto
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/scripts/kallisto/logs/dream_enhanced_%j.out
#SBATCH --error=RNA-seq/scripts/kallisto/logs/dream_enhanced_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/scripts/kallisto/logs

micromamba run -n rnaseq Rscript RNA-seq/scripts/kallisto/04b_dream_kallisto_enhanced.R
