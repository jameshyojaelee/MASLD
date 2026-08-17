#!/bin/bash
#SBATCH --job-name=fig3e_compact
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3e_compact_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3e_compact_%j.err

set -euo pipefail
: "${FIGURE_CANDIDATE_ROOT:?Set candidate root}"
: "${STAGE_RELEASE_ROOT:?Set validated stage release root}"
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${project_root}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
"${rnaseq_env}/bin/Rscript" --vanilla \
  "${project_root}/scripts/figures/progression_cascade.R"
pdf="${FIGURE_CANDIDATE_ROOT}/figure3/panels/fig3e_stage_remodeling.pdf"
test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
if strings "${pdf}" | grep -q '/Subtype /Type3'; then exit 1; fi
strings "${pdf}" | grep -q '/FontFamily (Helvetica)'
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
mkdir -p "${FIGURE_CANDIDATE_ROOT}/proofs"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r240 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/fig3e_stage_remodeling_compact.png" "${pdf}"
"${rnaseq_env}/bin/python3" \
  "${project_root}/scripts/figures/validate_fig3e_compact.py"
"${rnaseq_env}/bin/python3" \
  "${project_root}/scripts/manuscript/validate_resource_scope.py"
