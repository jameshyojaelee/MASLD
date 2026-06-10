#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/deconv_c2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/deconv_c2_%j.err

# Activate env BEFORE enabling nounset: conda/binutils activation scripts
# reference unbound vars (ADDR2LINE) and abort under `set -u`.
source ~/.bashrc
micromamba activate rnaseq
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "==================================================================="
echo "[C2 deconv recount] start $(date)"
echo "==================================================================="

echo ">>> [1/3] Script 25 C2: deconvolution attribution (hep-intrinsic) ..."
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/25_deconv_attribution_C2.R
echo ">>> [1/3] Script 25 C2 done $(date)"

echo ">>> [2/3] Script 80 C2: MuSiC signed-score cross-attribution ..."
Rscript RNA-seq/80_celltype_intrinsic_attribution_C2.R
echo ">>> [2/3] Script 80 C2 done $(date)"

echo ">>> [3/4] Script 89 C2: TOAST CT-specific DE (reproducibility) ..."
Rscript RNA-seq/89_toast_celltype_de_C2.R
echo ">>> [3/4] Script 89 C2 done $(date)"

echo ">>> [4/4] Script 80b: cross-method hepatocyte-intrinsic Jaccard (dream vs C2) ..."
Rscript RNA-seq/80b_c2_crossmethod_jaccard.R
echo ">>> [4/4] Script 80b done $(date)"

echo "==================================================================="
echo "[C2 deconv recount] ALL DONE $(date)"
echo "==================================================================="
