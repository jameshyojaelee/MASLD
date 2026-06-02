#!/bin/bash
#SBATCH --job-name=snapatac2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_Multiome/logs/corrected_da_%j.out
#SBATCH --error=Analysis/ATAC/Human_Multiome/logs/corrected_da_%j.err

# NB: no `-u` -- the rnaseq conda activate scripts reference unbound vars
# (e.g. ADDR2LINE) and abort Step 2 under `set -u`. (Fix 2026-05-30.)
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"
# NB: do NOT set PYTHONNOUSERSITE here. Step 1 (snapatac2) needs the ~/.local
# scipy for backed-anndata sparse indexing; hiding it breaks adata.X[idx,:].
# PYTHONNOUSERSITE=1 is set ONLY for Step 2 (rnaseq) below, where the ~/.local
# numexpr is broken. (Two steps need opposite site-package isolation.)
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2

ATAC_DIR="Analysis/ATAC/Human_Multiome"
RESULTS="${ATAC_DIR}/results/snapatac2"

echo "=== Step 1: Corrected DA (donor pseudobulk OLS, hepatocyte-only) ==="
python "${ATAC_DIR}/scripts/07b_corrected_da_hepatocytes.py" \
    --input "${ATAC_DIR}/results/label_transfer/snapatac2_label_transferred.h5ad" \
    --per-donor-dir "${RESULTS}/per_donor" \
    --output-dir "${RESULTS}" \
    --min-detection 0.05

echo ""
echo "=== Step 2: Annotate peaks with gene symbols ==="
# Ignore ~/.local here: its numexpr has a numpy-ABI break that crashes pandas
# import in rnaseq. (Scoped to Step 2 only; Step 1 needs ~/.local scipy.)
export PYTHONNOUSERSITE=1
micromamba activate rnaseq
python "${ATAC_DIR}/scripts/08_annotate_peaks_for_l8.py" \
    --da-results "${RESULTS}/scatac_da_corrected_hep.csv" \
    --chromvar "${ATAC_DIR}/results/chromvar/chromvar_tf_activity.csv" \
    --regulons "${ATAC_DIR}/scenic_plus/hepatocyte_regulons.csv" \
    --gtf /gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
    --output-dir "${ATAC_DIR}/results/l8_annotated_corrected"

echo ""
echo "=== Step 3 (figure regen) SKIPPED: figures are regenerated in the Track D"
echo "    wave from the corrected outputs; figS_atac.R may be renamed/stale. ==="

echo ""
echo "=== Done: $(date) ==="
ls -la "${RESULTS}/scatac_da_corrected_hep.csv" || true
ls -la "${ATAC_DIR}/results/l8_annotated_corrected/" || true
# (Step-3 volcano PDF removed; figures are regenerated in the Track D wave.)
