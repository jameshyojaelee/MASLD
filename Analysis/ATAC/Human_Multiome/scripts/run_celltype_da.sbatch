#!/bin/bash
#SBATCH --job-name=snapatac2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_Multiome/logs/celltype_da_%j.out
#SBATCH --error=Analysis/ATAC/Human_Multiome/logs/celltype_da_%j.err

# Donor-level pseudobulk DA for the 3 non-hepatocyte Fig S4g compartments, then
# edgeR-QLF, then regenerate Fig S4g from the 4-compartment donor-level set.
# NB: no `set -u` — rnaseq activate scripts reference unbound vars. (Same as
# run_corrected_da.sbatch.)
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$(pwd)"

eval "$(micromamba shell hook -s bash)"

ATAC_DIR="Analysis/ATAC/Human_Multiome"
RESULTS="${ATAC_DIR}/results/snapatac2"

echo "=== Step 1: donor pseudobulk (Stellate / Macrophage+Kupffer / Cholangiocyte) ==="
# snapatac2 env; needs ~/.local scipy for backed-anndata sparse indexing, so do
# NOT set PYTHONNOUSERSITE here.
micromamba activate snapatac2
python "${ATAC_DIR}/scripts/07c_corrected_da_celltypes.py" \
    --input "${ATAC_DIR}/results/label_transfer/snapatac2_label_transferred.h5ad" \
    --output-dir "${RESULTS}"

echo ""
echo "=== Step 2: edgeR-QLF per compartment (canonical donor-level test) ==="
# rnaseq env; ~/.local numexpr is ABI-broken -> hide it for the R/python steps.
export PYTHONNOUSERSITE=1
micromamba activate rnaseq
Rscript "${ATAC_DIR}/scripts/17c_edger_celltypes.R"

echo ""
echo "=== Step 3: regenerate Fig S4g from the 4-compartment donor-level set ==="
Rscript scripts/figures/figS05_scatac_disease_vs_control.R

echo ""
echo "=== Done: $(date) ==="
ls -la "${RESULTS}"/scatac_da_corrected_*_edger.csv || true
ls -la figures/main/fig4_validation/panels/figS4g.pdf || true
