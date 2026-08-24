#!/usr/bin/env bash
#SBATCH --job-name=hac_tx_plot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_ANALYSIS_CANDIDATE:?HAC_ANALYSIS_CANDIDATE is unset}"
: "${HAC_FIGURE_CANDIDATE:?HAC_FIGURE_CANDIDATE is unset}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT_DIR}/61_render_transcript_overview.R')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT_DIR}/61_render_transcript_overview.R" \
  "${HAC_ANALYSIS_CANDIDATE}" "${HAC_FIGURE_CANDIDATE}"
pdfinfo "${HAC_FIGURE_CANDIDATE}/panels/s3_continuum_transcript_coprimary_agreement.pdf"
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage \
  "${HAC_FIGURE_CANDIDATE}/panels/s3_continuum_transcript_coprimary_agreement.pdf"
sha256sum "${HAC_FIGURE_CANDIDATE}/panels/"*.pdf > "${HAC_FIGURE_CANDIDATE}/checksums.sha256"
