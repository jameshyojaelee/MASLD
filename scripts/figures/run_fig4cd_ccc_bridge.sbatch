#!/bin/bash
#SBATCH --job-name=fig4cd_ccc_bridge
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4cd_ccc_bridge_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4cd_ccc_bridge_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE_ROOT="${PROJECT_ROOT}/figures/candidates/fig4cd-ccc-bridge-2026-08-15-v10"
PROGRAM_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"
CURRENT_4C="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/panels/fig4c_hotspot_stage_heatmap.pdf"
RNASEQ_RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

test -s "${PROGRAM_ROOT}/READY"
test -s "${CURRENT_4C}"
test ! -e "${CANDIDATE_ROOT}"
mkdir -p "${CANDIDATE_ROOT}/proofs"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export FIGURE_CANDIDATE_ROOT="${CANDIDATE_ROOT}"
export PROGRAM_RELEASE_ROOT="${PROGRAM_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RNASEQ_RSCRIPT}" --vanilla scripts/figures/render_fig4cd_ccc_bridge.R
"${RNASEQ_RSCRIPT}" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIGURE_CANDIDATE_ROOT"), "provenance/sessionInfo.txt"))'

python3 scripts/manuscript/validate_resource_scope.py

install -m 0444 scripts/figures/render_fig4cd_ccc_bridge.R \
  "${CANDIDATE_ROOT}/provenance/render_fig4cd_ccc_bridge.R"
install -m 0444 scripts/figures/validate_fig4cd_ccc_bridge.py \
  "${CANDIDATE_ROOT}/provenance/validate_fig4cd_ccc_bridge.py"
install -m 0444 scripts/figures/run_fig4cd_ccc_bridge.sbatch \
  "${CANDIDATE_ROOT}/provenance/run_fig4cd_ccc_bridge.sbatch"

for pdf in "${CANDIDATE_ROOT}"/panels/*.pdf; do
  panel="$(basename "${pdf}" .pdf)"
  gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r180 \
    -sOutputFile="${CANDIDATE_ROOT}/proofs/${panel}.png" "${pdf}"
done
montage \
  "${CANDIDATE_ROOT}/proofs/fig4c_hotspot_stage_heatmap.png" \
  "${CANDIDATE_ROOT}/proofs/fig4d_program_linked_communication.png" \
  -tile 2x1 -geometry +24+24 -background white \
  "${CANDIDATE_ROOT}/proofs/fig4cd_bridge_contact_sheet.png"

python3 scripts/figures/validate_fig4cd_ccc_bridge.py \
  "${CANDIDATE_ROOT}" "${CURRENT_4C}"
