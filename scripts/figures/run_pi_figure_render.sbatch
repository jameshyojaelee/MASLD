#!/bin/bash
#SBATCH --job-name=pi_fig_render
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/pi_fig_render_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/pi_fig_render_%j.err

set -euo pipefail
umask 027

: "${FIGURE_CANDIDATE_ROOT:?Set the candidate path created by the analysis job}"
: "${STAGE_RELEASE_ROOT:?Set the validated five-cohort fragment release}"
test -d "${FIGURE_CANDIDATE_ROOT}/analysis/stage_extensions"
test ! -e "${FIGURE_CANDIDATE_ROOT}/output_manifest.tsv"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export PATH="${RNASEQ_ENV}/bin:/gpfs/commons/home/jameslee/.local/bin:${PATH}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

R="${RNASEQ_ENV}/bin/Rscript"
"${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/nas_fib_grid.R"
"${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/figS_pca_definitive.R"
"${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/progression_cascade.R"
"${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/gen_scrna_umap_embeddable.R"
"${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/singlecell_module_heatmap.R"
"${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/tf_convergence_vsF0.R"
"${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/fig5a_input_firewall.R"
"${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/fig5b_protein_triage.R"
"${R}" --vanilla -e 'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIGURE_CANDIDATE_ROOT"), "sessionInfo.txt"))'

if [[ -n "${GENETICS_RELEASE_ROOT:-}" ]]; then
  "${R}" --vanilla "${PROJECT_ROOT}/scripts/figures/stage_deg_genetics_matrix.R"
fi

"${RNASEQ_ENV}/bin/python3" "${PROJECT_ROOT}/scripts/figures/finalize_pi_figure_candidate.py"
"${RNASEQ_ENV}/bin/python3" "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
