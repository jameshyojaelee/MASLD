#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/primary_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/primary_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2
ENV=${ROOT}/.mamba/seqfunc_leafcutter
export PYTHONNOUSERSITE=1
micromamba run -p "${ENV}" leafcutter-ds "${OUT}/joint/joint846_perind_numers.counts.gz" "${OUT}/joint/groups.tsv" \
  -0 Control -o "${OUT}/joint/primary_dataset_sex" -i 5 -g 10 -c 20 -p 32 -e "${OUT}/gencode_v49_exons.tsv.gz"
test -s "${OUT}/joint/primary_dataset_sex_cluster_significance.txt"
test -s "${OUT}/joint/primary_dataset_sex_effect_sizes.txt"
