#!/usr/bin/env bash
#SBATCH --job-name=fig2panel
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/fig2panel_%j.out
#SBATCH --error=logs/fig2panel_%j.err

set -euo pipefail
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rscript_bin="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq_r44/bin/Rscript"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"

: "${FIG2_TRAIT_DIRECTNESS_STATS:?set FIG2_TRAIT_DIRECTNESS_STATS to the stats candidate dir}"
export MASLD_PROJECT_ROOT="${project_root}"
export FIG2_TRAIT_DIRECTNESS_OUT="${project_root}/figures/candidates/fig2_trait_directness_panels_${stamp}_${SLURM_JOB_ID}"
mkdir -p "${project_root}/logs"
cd "${project_root}"

test ! -e "${FIG2_TRAIT_DIRECTNESS_OUT}"
mkdir -p "${FIG2_TRAIT_DIRECTNESS_OUT}"

"${rscript_bin}" scripts/figures/fig2_trait_directness_panels.R \
  > "${FIG2_TRAIT_DIRECTNESS_OUT}/render.log" 2>&1

for pdf in "${FIG2_TRAIT_DIRECTNESS_OUT}"/*.pdf; do
  name="$(basename "${pdf}" .pdf)"
  pdfinfo "${pdf}" > "${FIG2_TRAIT_DIRECTNESS_OUT}/${name}_pdfinfo.txt"
  pdftotext "${pdf}" "${FIG2_TRAIT_DIRECTNESS_OUT}/${name}_artwork_text.txt"
  pages="$(awk '/^Pages:/ {print $2}' "${FIG2_TRAIT_DIRECTNESS_OUT}/${name}_pdfinfo.txt")"
  if [ "${pages}" != "1" ]; then
    echo "expected one page in ${name}, got ${pages}" >&2
    exit 1
  fi
  if grep -Eqi 'causal gene|causal variant|proves|longitudinal progression' \
    "${FIG2_TRAIT_DIRECTNESS_OUT}/${name}_artwork_text.txt"; then
    echo "prohibited wording in ${name}" >&2
    exit 1
  fi
done

sha256sum "${FIG2_TRAIT_DIRECTNESS_OUT}"/*.pdf "${FIG2_TRAIT_DIRECTNESS_OUT}"/*.tsv \
  > "${FIG2_TRAIT_DIRECTNESS_OUT}/checksums.txt"
printf '%s\n' "${FIG2_TRAIT_DIRECTNESS_OUT}"
