#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --array=0-9%5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/perm_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/perm_%A_%a.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2
ENV=${ROOT}/.mamba/seqfunc_leafcutter
i=${SLURM_ARRAY_TASK_ID}; tag=$(printf '%02d' "$i")
export PYTHONNOUSERSITE=1
micromamba run -p "${ENV}" leafcutter-ds "${OUT}/joint/joint846_perind_numers.counts.gz" "${OUT}/permutations/groups_perm_${tag}.tsv" \
  -0 Control -o "${OUT}/permutations/perm_${tag}" -i 5 -g 10 -c 20 -p 24
test -s "${OUT}/permutations/perm_${tag}_cluster_significance.txt"
