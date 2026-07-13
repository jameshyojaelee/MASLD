#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_joint_ds_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_joint_ds_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_leafcutter
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/production/joint_native
export PYTHONNOUSERSITE=1
counts=${OUT}/joint860_perind_numers.counts.gz
test -s "${counts}"
micromamba run -p "${ENV}" leafcutter-ds "${counts}" "${OUT}/groups_dataset_fixed.tsv" \
  -0 Control -o "${OUT}/joint860_dataset_fixed" -i 5 -g 10 -c 20 -p 16 \
  -e "${OUT}/gencode_v49_exons.tsv.gz"
test -s "${OUT}/joint860_dataset_fixed_cluster_significance.txt"
test -s "${OUT}/joint860_dataset_fixed_effect_sizes.txt"
