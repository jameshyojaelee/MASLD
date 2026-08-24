#!/bin/bash
#SBATCH --job-name=fig3g_promoted
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3g_promoted_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3g_promoted_%j.err

set -euo pipefail
umask 027
: "${FIGURE_CANDIDATE_ROOT:?Set a fresh candidate path}"
: "${STAGE_RELEASE_ROOT:?Set the validated fragment release}"
: "${GENETICS_RELEASE_ROOT:?Set the promoted corrected-COLOC release}"
: "${GENETICS_STATE_INTERFACE_ROOT:?Set the promoted-release interface result}"
test ! -e "${FIGURE_CANDIDATE_ROOT}"
mkdir -p "${FIGURE_CANDIDATE_ROOT}"
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${project_root}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
"${rnaseq_env}/bin/Rscript" --vanilla \
  "${project_root}/scripts/figures/fig3g_genetic_state_complementarity.R"
pdf="${FIGURE_CANDIDATE_ROOT}/figure3/panels/fig3g_genetic_state_complementarity.pdf"
supp_pdf="${FIGURE_CANDIDATE_ROOT}/figureS3/panels/figs3g_expression_matched_pairs.pdf"
test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
test "$(pdfinfo "${supp_pdf}" | awk '/^Pages:/ {print $2}')" = "1"
if strings "${pdf}" | grep -q '/Subtype /Type3'; then exit 1; fi
if strings "${supp_pdf}" | grep -q '/Subtype /Type3'; then exit 1; fi
strings "${pdf}" | grep -q '/FontFamily (Helvetica)'
strings "${supp_pdf}" | grep -q '/FontFamily (Helvetica)'
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${supp_pdf}"
awk -F '\t' 'NR == 2 {if ($2 != 14112 || $3 != 3617 || $7 != 428 ||
  $8 != 30 || $9 != 428 || $10 != 30 || $11 != 21 || $12 != 377 ||
  $13 != 0 || $16 != 0.25) exit 1}' \
  "${FIGURE_CANDIDATE_ROOT}/source_tables/fig3g_matched_treat_summary.tsv"
mkdir -p "${FIGURE_CANDIDATE_ROOT}/proofs"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r240 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/fig3g_genetic_state_complementarity.png" "${pdf}"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r240 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/figs3g_expression_matched_pairs.png" "${supp_pdf}"
find "${FIGURE_CANDIDATE_ROOT}" -type f ! -name input_output_manifest.tsv -print0 | \
  sort -z | xargs -0 sha256sum > "${FIGURE_CANDIDATE_ROOT}/input_output_manifest.tsv"
"${rnaseq_env}/bin/python3" \
  "${project_root}/scripts/manuscript/validate_resource_scope.py"
echo "FIG3G_PROMOTED_RELEASE_GATE PASS"
