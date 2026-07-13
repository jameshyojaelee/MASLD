#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=72:00:00
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/cohort_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/cohort_%A_%a.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2
ENV=${ROOT}/.mamba/seqfunc_leafcutter
cohorts=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)
ds=${cohorts[$SLURM_ARRAY_TASK_ID]}; d=${OUT}/cohorts/${ds}
export PYTHONNOUSERSITE=1
micromamba run -p "${ENV}" leafcutter-cluster -j "${d}/juncfiles.txt" -r "${d}" -o "${ds}" -l 500000 -m 30 -p 0.001
micromamba run -p "${ENV}" leafcutter-ds "${d}/${ds}_perind_numers.counts.gz" "${d}/groups.tsv" \
  -0 Control -o "${d}/${ds}_sex" -i 5 -g 10 -c 20 -p 16 -e "${OUT}/gencode_v49_exons.tsv.gz"
test -s "${d}/${ds}_sex_cluster_significance.txt"
