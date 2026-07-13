#!/bin/bash
#SBATCH --job-name=diffvar
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --output=scripts/heterogeneity_program/phase0/logs/t1_dv_%j.out
set -eo pipefail   # NOT -u: rnaseq binutils activate script has unbound ADDR2LINE
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/heterogeneity_program/phase0/logs
source ~/.bashrc
micromamba activate rnaseq
export N_PERM=1000
Rscript scripts/heterogeneity_program/phase0/t1_dv_screen.R
