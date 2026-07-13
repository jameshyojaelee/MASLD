#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/loco_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/loco_%A_%a.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2
ENV=${ROOT}/.mamba/seqfunc_leafcutter
cohorts=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)
ds=${cohorts[$SLURM_ARRAY_TASK_ID]}; d=${OUT}/loco/${ds}
export PYTHONNOUSERSITE=1
micromamba run -p "${ENV}" leafcutter-ds "${OUT}/joint/joint846_perind_numers.counts.gz" "${d}/groups.tsv" \
  -0 Control -o "${d}/loco_${ds}" -i 5 -g 10 -c 20 -p 24 -e "${OUT}/gencode_v49_exons.tsv.gz"
test -s "${d}/loco_${ds}_cluster_significance.txt"
