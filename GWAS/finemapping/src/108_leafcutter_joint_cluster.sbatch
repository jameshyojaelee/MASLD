#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_joint_cluster_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_joint_cluster_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_leafcutter
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/production/joint_native
export PYTHONNOUSERSITE=1 MASLD_PROJECT_ROOT=${ROOT}
python3 "${ROOT}/GWAS/finemapping/src/108_leafcutter_joint_prepare.py"
test "$(wc -l < "${OUT}/juncfiles.txt")" -eq 860
micromamba run -p "${ENV}" leafcutter-cluster -j "${OUT}/juncfiles.txt" -r "${OUT}" -o joint860 -l 500000 -m 30 -p 0.001
test -s "${OUT}/joint860_perind_numers.counts.gz"
test -s "${OUT}/joint860_perind.counts.gz"
