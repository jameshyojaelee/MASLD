#!/usr/bin/env bash
#SBATCH --job-name=hac_competitive_null
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/74_render_system_competitive_null.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_SYSTEM_STRESS_SOURCE:?HAC_SYSTEM_STRESS_SOURCE is unset}"
: "${HAC_COMPETITIVE_NULL_CANDIDATE:?HAC_COMPETITIVE_NULL_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" \
  "${HAC_SYSTEM_STRESS_SOURCE}" "${HAC_COMPETITIVE_NULL_CANDIDATE}"

for PANEL in \
  "${HAC_COMPETITIVE_NULL_CANDIDATE}/panels/s4v_system01_competitive_null.pdf" \
  "${HAC_COMPETITIVE_NULL_CANDIDATE}/panels/s4w_system_competitive_gate.pdf"; do
  [[ "$(pdfinfo "${PANEL}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PANEL}"
  ! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${PANEL}"
done

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_COMPETITIVE_NULL_CANDIDATE")
  m <- fread(file.path(root, "figure_manifest.tsv"))
  g <- fread(file.path(root, "source_tables", "s4w_system_competitive_gate.tsv"))
  f <- fread(file.path(root, "source_tables", "s4v_system01_competitive_null.tsv"))
  a <- fread(file.path(root, "provenance", "audit.tsv"))
  stopifnot(
    nrow(m) == 2L, all(m$callout == c("S4V", "S4W")),
    all(m$connected_null_supported == "0/43"),
    all(m$uniform_null_supported == "1/43"),
    nrow(g) == 86L, nrow(f) == 4L,
    sum(g[null_model == "connected_subgraph", gate_q_value < 0.05]) == 0L,
    sum(g[null_model == "size_matched_uniform", gate_q_value < 0.05]) == 1L,
    all(a$passed)
  )
  cat("competitive-null panel assertions passed\n")
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
