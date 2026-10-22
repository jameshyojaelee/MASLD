#!/bin/bash
#SBATCH --job-name=9522
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F2abce/jobs/9522/%j.out
for f in /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/combined_finemapping.csv /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/combined_finemapping_coloc_targeted.csv; do echo "== $f"; grep -E '^1,155086571,' $f | grep UKBB_GGT; grep -E '^1,155129562,' $f | grep UKBB_ALT; done
