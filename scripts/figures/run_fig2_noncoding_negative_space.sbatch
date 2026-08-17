#!/usr/bin/env bash
#SBATCH --job-name=fig2_negspace
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=logs/fig2_negspace_%j.out
#SBATCH --error=logs/fig2_negspace_%j.err

set -euo pipefail
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rscript_bin="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq_r44/bin/Rscript"
stamp="$(date +%Y%m%d_%H%M%S)"
export MASLD_PROJECT_ROOT="${project_root}"
export FIG2_NEGATIVE_SPACE_DIR="${project_root}/figures/candidates/fig2_noncoding_negative_space_${stamp}_${SLURM_JOB_ID}"
mkdir -p "${project_root}/logs" "${FIG2_NEGATIVE_SPACE_DIR}"
cd "${project_root}"

"${rscript_bin}" scripts/figures/fig2_noncoding_negative_space.R \
  > "${FIG2_NEGATIVE_SPACE_DIR}/render.log" 2>&1

pdf="${FIG2_NEGATIVE_SPACE_DIR}/Fig2_companion_noncoding_negative_space.pdf"
pdfinfo "${pdf}" > "${FIG2_NEGATIVE_SPACE_DIR}/pdfinfo.txt"
pdftotext "${pdf}" "${FIG2_NEGATIVE_SPACE_DIR}/artwork_text.txt"
if grep -Eqi 'causal gene|ancestry-specific|ABF-fallback|exploratory' \
  "${FIG2_NEGATIVE_SPACE_DIR}/artwork_text.txt"; then
  echo "prohibited wording in artwork" >&2
  exit 1
fi
sha256sum \
  figures/candidates/noncoding-resource-2026-08-12/noncoding-dna-v5/fig2_noncoding_regulatory_context_source.tsv \
  > "${FIG2_NEGATIVE_SPACE_DIR}/input_sha256.txt"
printf '%s\n' "${FIG2_NEGATIVE_SPACE_DIR}"
