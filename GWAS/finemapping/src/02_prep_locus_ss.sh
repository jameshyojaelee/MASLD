#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu

# Prep per-locus summary stats SLURM wrapper
# Usage: sbatch 02_prep_locus_ss.sh <sumstats_name> <sumstats_path> <leadsnps_path> <ld_pop> <window_mb>

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping

eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping

Rscript src/02_prep_locus_ss.R "$1" "$2" "$3" "$4" "$5"
