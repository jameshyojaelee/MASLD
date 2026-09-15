#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_atac_allelic
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/logs/atac-allelic-%j.out
set -euo pipefail
bash /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas/72_atac_allelic_counts.sh /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_External/cellranger /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/atac_counts /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/atac_regions.txt 12
