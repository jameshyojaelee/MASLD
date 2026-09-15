#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_atac_allelic_244832
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/logs/atac-allelic-244832-%j.out
set -euo pipefail
P=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
bash $P/scripts/analysis/alphagenome_atlas/72_atac_allelic_counts.sh \
  $P/Analysis/ATAC/Human_Multiome/results/alignment \
  $P/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/atac_counts_gse244832 \
  $P/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/atac_regions.txt \
  16 flat
