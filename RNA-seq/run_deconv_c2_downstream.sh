#!/bin/bash
#SBATCH --job-name=MuSiC
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/deconv_c2_down_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/deconv_c2_down_%j.err

# C2 deconv recount — DOWNSTREAM ONLY.
# Script 25_C2 already wrote the primary MuSiC attribution
# (deconv_attribution_scores.csv: Hep_intrinsic=1135 / Composition=718 on the
# 1,853 C2 Tier-1) before crashing in the optional legacy 2-CT side table.
# These three steps consume the primary attribution + scRNA + bulk DGE and do
# NOT need the expensive dream fits redone.

# Activate env BEFORE nounset (conda binutils activation references unbound vars).
source ~/.bashrc
micromamba activate rnaseq
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "[C2 downstream] start $(date)"

echo ">>> [1/3] Script 80 C2: MuSiC signed-score cross-attribution ..."
Rscript RNA-seq/80_celltype_intrinsic_attribution_C2.R
echo ">>> [1/3] Script 80 C2 done $(date)"

echo ">>> [2/3] Script 89 C2: TOAST CT-specific DE (reproducibility) ..."
Rscript RNA-seq/89_toast_celltype_de_C2.R
echo ">>> [2/3] Script 89 C2 done $(date)"

echo ">>> [3/3] Script 80b: cross-method hepatocyte-intrinsic Jaccard (dream vs C2) ..."
Rscript RNA-seq/80b_c2_crossmethod_jaccard.R
echo ">>> [3/3] Script 80b done $(date)"

echo "[C2 downstream] ALL DONE $(date)"
