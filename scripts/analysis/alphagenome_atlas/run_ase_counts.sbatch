#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_ase_GSE126848
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/logs/ase-GSE126848-%j.out
set -euo pipefail
bash /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas/70_ase_counts.sh GSE126848 /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE126848/alignments/star /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/ase_targets.tsv 16
