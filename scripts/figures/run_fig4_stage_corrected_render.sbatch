#!/bin/bash
#SBATCH --job-name=fig4_stage_render
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_stage_render_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_stage_render_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FIGURE_ROOT="${PROJECT_ROOT}/figures/candidates/fig4-stage-terminology-corrected-2026-08-13-v4"
PROGRAM_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3"
STAGE_ROOT="${PROJECT_ROOT}/figures/candidates/pi-figure-redesign-2026-08-13-v6/analysis/stage_extensions"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"

test -s "${PROGRAM_ROOT}/hotspot/READY"
test -s "${STAGE_ROOT}/output_manifest.tsv"
test -s "${FIGURE_ROOT}/source_tables/fig4d_communication_trajectories.tsv"
test -s "${FIGURE_ROOT}/source_tables/fig4d_communication_audit.tsv"
test ! -e "${FIGURE_ROOT}/figure4/panels/fig4b_scrna_umap_embeddable.pdf"
test ! -e "${FIGURE_ROOT}/figure4/panels/fig4c_hotspot_stage_heatmap.pdf"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export FIGURE_CANDIDATE_ROOT="${FIGURE_ROOT}"
export STAGE_RELEASE_ROOT="${STAGE_ROOT}"
export PROGRAM_RELEASE_ROOT="${PROGRAM_ROOT}/hotspot"
export SINGLECELL_STAGE_METADATA_FILE="${PROGRAM_ROOT}/metadata/donor_metadata_stage_corrected.tsv"
export CCC_TRAJECTORY_FILE="${FIGURE_ROOT}/source_tables/fig4d_communication_trajectories.tsv"
export CCC_AUDIT_FILE="${FIGURE_ROOT}/source_tables/fig4d_communication_audit.tsv"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RNASEQ_ENV}/bin/Rscript" --vanilla scripts/figures/gen_scrna_umap_embeddable.R
"${RNASEQ_ENV}/bin/Rscript" --vanilla scripts/figures/singlecell_module_heatmap.R
"${RNASEQ_ENV}/bin/Rscript" --vanilla scripts/figures/render_pi_singlecell_panels.R
"${RNASEQ_ENV}/bin/Rscript" --vanilla scripts/figures/tf_convergence_vsF0.R
"${RNASEQ_ENV}/bin/Rscript" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIGURE_CANDIDATE_ROOT"), "sessionInfo.txt"))'

python3 scripts/manuscript/validate_resource_scope.py
