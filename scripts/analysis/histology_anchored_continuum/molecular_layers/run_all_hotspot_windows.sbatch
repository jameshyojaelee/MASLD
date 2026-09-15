#!/usr/bin/env bash
#SBATCH --job-name=hac_hotspot_windows
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_ANALYSIS_CANDIDATE:?HAC_ANALYSIS_CANDIDATE is unset}"
: "${HAC_HOTSPOT_WINDOW_CANDIDATE:?HAC_HOTSPOT_WINDOW_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT_DIR}/63_render_all_hotspot_windows.R')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT_DIR}/63_render_all_hotspot_windows.R" \
  "${HAC_ANALYSIS_CANDIDATE}" "${HAC_HOTSPOT_WINDOW_CANDIDATE}"

while IFS=$'\t' read -r artifact path _; do
  [[ "${artifact}" == "artifact" ]] && continue
  [[ "$(pdfinfo "${path}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${path}"
  if grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${path}"; then
    echo "Type 3 font detected: ${path}" >&2
    exit 1
  fi
done < "${HAC_HOTSPOT_WINDOW_CANDIDATE}/figure_manifest.tsv"

sha256sum -c <(
  awk 'BEGIN{FS="\t"} NR>1 {print $8 "  " $2}' \
    "${HAC_HOTSPOT_WINDOW_CANDIDATE}/figure_manifest.tsv"
)

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_HOTSPOT_WINDOW_CANDIDATE")
  x <- fread(file.path(root, "source_tables", "all117_hotspot_continuum_windows.tsv.gz"))
  s <- fread(file.path(root, "summary.tsv"))
  stopifnot(uniqueN(x$program_uid) == 117L, uniqueN(x$dataset) == 2L,
            uniqueN(x$window_id) == 9L, nrow(x) == 117L * 2L * 9L,
            all(x$visualization_only),
            s[metric == "continuum_associated_programs", as.integer(value)] > 0L)
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
