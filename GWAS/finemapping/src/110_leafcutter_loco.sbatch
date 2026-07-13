#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_loco_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_loco_%A_%a.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_leafcutter
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/production/joint_native
cohorts=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)
drop=${cohorts[$SLURM_ARRAY_TASK_ID]}
export PYTHONNOUSERSITE=1
micromamba run -p "${ENV}" leafcutter-ds "${OUT}/joint860_perind_numers.counts.gz" \
  "${OUT}/groups_loco_${drop}.tsv" -0 Control -o "${OUT}/loco_${drop}" \
  -i 5 -g 10 -c 20 -p 16 -e "${OUT}/gencode_v49_exons.tsv.gz"
test -s "${OUT}/loco_${drop}_cluster_significance.txt"
