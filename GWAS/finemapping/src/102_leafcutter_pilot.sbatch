#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_leafcutter
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/pilot
export PYTHONNOUSERSITE=1
test "$(find "${OUT}/junctions" -name '*.junc' -size +0c | wc -l)" -eq 20
counts=${OUT}/pilot_perind_numers.counts.gz
if [ ! -s "${counts}" ]; then
  micromamba run -p "${ENV}" leafcutter-cluster -j "${OUT}/juncfiles.txt" -r "${OUT}" -o pilot -l 500000 -m 30 -p 0.001
fi
test -s "${counts}"
micromamba run -p "${ENV}" leafcutter-ds "${counts}" "${OUT}/groups.tsv" -0 Control -o "${OUT}/pilot_ds" -i 3 -g 2 -c 20 -p 8
test -s "${OUT}/pilot_ds_cluster_significance.txt"
