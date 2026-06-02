#!/usr/bin/env bash
#SBATCH --job-name=export_cibersortx
#SBATCH --cpus-per-task=8
#SBATCH --mem=600G
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem
#SBATCH --output=RNA-seq/deconvolution/logs/cibersortx_export_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/cibersortx_export_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13

# Environment setup
export LD_LIBRARY_PATH="RNA-seq/.mamba/music_deconv/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export R_LIBS_USER=
export R_LIBS=
export R_LIBS_SITE="RNA-seq/.mamba/music_deconv/lib/R/library"

MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv

REF_MOUSE=${REF_MOUSE:-RNA-seq/deconvolution/reference/reference_mouse_sce.rds}
REF_HUMAN=${REF_HUMAN:-RNA-seq/deconvolution/reference/reference_human_sce.rds}

OUTPUT_DIR=RNA-seq/deconvolution/cibersortx_input
mkdir -p "${OUTPUT_DIR}"

echo "=== Exporting Mouse reference for CIBERSORTx ==="
"${MICROMAMBA}" run -p "${ENV_PREFIX}" \
  Rscript --vanilla \
  RNA-seq/deconvolution/scripts/11_export_cibersortx.R \
  "${REF_MOUSE}" \
  "${OUTPUT_DIR}/mouse" \
  "celltype"

echo ""
echo "=== Exporting Human reference for CIBERSORTx ==="
"${MICROMAMBA}" run -p "${ENV_PREFIX}" \
  Rscript --vanilla \
  RNA-seq/deconvolution/scripts/11_export_cibersortx.R \
  "${REF_HUMAN}" \
  "${OUTPUT_DIR}/human" \
  "celltype"

echo ""
echo "=== CIBERSORTx Export Complete ==="
echo "Output files in: ${OUTPUT_DIR}"
echo ""
echo "Next steps:"
echo "1. Go to https://cibersortx.stanford.edu/"
echo "2. Upload *_sc_matrix.txt and *_phenoclasses.txt"
echo "3. Create signature matrix"
echo "4. Download signature matrix and use for deconvolution"
