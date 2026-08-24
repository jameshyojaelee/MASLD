#!/usr/bin/env bash
#SBATCH --job-name=hac_stage_compare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_ANALYSIS_CANDIDATE:?HAC_ANALYSIS_CANDIDATE is unset}"
: "${HAC_FIGURE_CANDIDATE:?HAC_FIGURE_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT_DIR}/62_render_stage_continuum_comparison.R')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT_DIR}/62_render_stage_continuum_comparison.R" \
  "${HAC_ANALYSIS_CANDIDATE}" "${HAC_FIGURE_CANDIDATE}"

while IFS=$'\t' read -r artifact path _; do
  [[ "${artifact}" == "artifact" ]] && continue
  pdfinfo "${path}" >/dev/null
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${path}"
done < "${HAC_FIGURE_CANDIDATE}/figure_manifest.tsv"

sha256sum -c <(
  awk 'BEGIN{FS="\t"} NR>1 {print $3 "  " $2}' \
    "${HAC_FIGURE_CANDIDATE}/figure_manifest.tsv"
)
