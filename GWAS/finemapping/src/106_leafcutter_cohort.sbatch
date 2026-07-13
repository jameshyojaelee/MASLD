#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_cohort_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_cohort_%j.err
set -euo pipefail
: "${DATASET:?export DATASET}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_leafcutter
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/production/${DATASET}
export PYTHONNOUSERSITE=1
expected=$(wc -l < "${OUT}/juncfiles.txt"); observed=$(find "${OUT}/junctions" -name '*.junc' -size +0c | wc -l)
test "${observed}" -eq "${expected}"
micromamba run -p "${ENV}" leafcutter-cluster -j "${OUT}/juncfiles.txt" -r "${OUT}" -o "${DATASET}" -l 500000 -m 30 -p 0.001
counts=${OUT}/${DATASET}_perind_numers.counts.gz; test -s "${counts}"
micromamba run -p "${ENV}" leafcutter-ds "${counts}" "${OUT}/groups.tsv" -0 Control -o "${OUT}/${DATASET}_ds" -i 5 -g 10 -c 20 -p 8
test -s "${OUT}/${DATASET}_ds_cluster_significance.txt"
