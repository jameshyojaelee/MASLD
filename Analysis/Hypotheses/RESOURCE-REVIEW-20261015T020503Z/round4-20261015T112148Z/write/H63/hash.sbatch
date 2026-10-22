#!/bin/bash
#SBATCH --job-name=8033
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H63/8033_%j.out
Y=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026
sha256sum $Y/v.mat $Y/zon_struct_all_full.mat $Y/human_samples_metadata.xlsx
stat -c '%s %n' $Y/v.mat
