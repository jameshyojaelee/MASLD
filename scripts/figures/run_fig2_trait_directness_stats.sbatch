#!/usr/bin/env bash
#SBATCH --job-name=fig2arch
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/fig2arch_%j.out
#SBATCH --error=logs/fig2arch_%j.err

set -euo pipefail
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
python_bin="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
out_dir="${project_root}/figures/candidates/fig2_trait_directness_${stamp}_${SLURM_JOB_ID}"

export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "${project_root}/logs"
cd "${project_root}"

test ! -e "${out_dir}"
"${python_bin}" -B scripts/figures/fig2_trait_directness_stats.py \
  --project-root "${project_root}" \
  --output-dir "${out_dir}" 2>&1 | tee "${project_root}/logs/fig2arch_stats_${SLURM_JOB_ID}.txt"

cp "${project_root}/logs/fig2arch_stats_${SLURM_JOB_ID}.txt" "${out_dir}/stats.log"
printf '%s\n' "${out_dir}"
