#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --array=0-859%40
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_convert_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/leafcutter_convert_%A_%a.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
MAN=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/production/production_manifest.tsv
line=$(awk -F '\t' -v i="${SLURM_ARRAY_TASK_ID}" 'NR>1 && $1==i {print;exit}' "${MAN}"); test -n "${line}"
src=$(cut -f5 <<<"${line}"); dst=$(cut -f6 <<<"${line}")
python3 "${ROOT}/GWAS/finemapping/src/105_star_sj_to_leafcutter.py" "${src}" "${dst}"
test -s "${dst}"
