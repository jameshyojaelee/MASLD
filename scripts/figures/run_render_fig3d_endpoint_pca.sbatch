#!/bin/bash
#SBATCH --job-name=fig3d_label
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3d_label_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3d_label_%j.err

set -euo pipefail
: "${FIGURE_CANDIDATE_ROOT:?Set candidate root}"
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export FIG3D_PC1_PVE=7.6
export FIG3D_PC2_PVE=3.9
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
"${rnaseq_env}/bin/Rscript" --vanilla \
  "${MASLD_PROJECT_ROOT}/scripts/figures/render_fig3d_endpoint_pca.R"
pdf="${FIGURE_CANDIDATE_ROOT}/figure3/panels/fig3d_pca_fibrosis_gradient.pdf"
test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
if strings "${pdf}" | grep -q '/Subtype /Type3'; then exit 1; fi
strings "${pdf}" | grep -q '/FontFamily (Helvetica)'
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
mkdir -p "${FIGURE_CANDIDATE_ROOT}/proofs"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r200 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/fig3d_endpoint_pca.png" "${pdf}"
