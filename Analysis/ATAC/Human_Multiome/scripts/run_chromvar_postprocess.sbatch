#!/bin/bash
#SBATCH --job-name=chromVAR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/chromvar_post_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/chromvar_post_%j.err

# =============================================================================
# Post-chromVAR: Re-run annotation (08) + L8 integration (35) + fig2
# Depends on: chromVAR re-run (03_chromvar_motifs.py)
# =============================================================================

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATAC_DIR="${PROJECT_ROOT}/Analysis/ATAC/Human_Multiome"

echo "============================================================"
echo "chromVAR Post-Processing Pipeline"
echo "  Job ID:     ${SLURM_JOB_ID}"
echo "  Started:    $(date)"
echo "============================================================"

module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

# Step 1: Verify chromVAR output exists and is corrected
CHROMVAR_CSV="${ATAC_DIR}/results/chromvar/chromvar_tf_activity.csv"
if [[ ! -f "$CHROMVAR_CSV" ]]; then
    echo "ERROR: chromVAR output not found: $CHROMVAR_CSV"
    exit 1
fi

# Quick sanity check: no corrupted motifs
n_extreme=$(python3 -c "
import pandas as pd
cv = pd.read_csv('$CHROMVAR_CSV')
hep = cv[cv.cell_type == 'Hepatocyte']
n_bad = ((hep.mean_deviation_masld.abs() > 5) | (hep.mean_deviation_normal.abs() > 5)).sum()
print(n_bad)
")
echo "Corrupted motifs in hepatocytes: $n_extreme"
if [[ "$n_extreme" -gt 0 ]]; then
    echo "WARNING: Still have $n_extreme corrupted motifs. Check chromVAR run."
fi

# Step 2: Re-run peak annotation for L8
echo ""
echo "Step 2: Re-annotating chromVAR for L8..."
python "${ATAC_DIR}/scripts/08_annotate_peaks_for_l8.py" \
    --da-results "${ATAC_DIR}/results/snapatac2/scatac_da_results.csv" \
    --chromvar "${CHROMVAR_CSV}" \
    --regulons "${ATAC_DIR}/scenic_plus/hepatocyte_regulons.csv" \
    --gtf "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz" \
    --output-dir "${ATAC_DIR}/results/l8_annotated"
echo "  Annotation done"

# Step 3: Re-run L8 ATAC integration into atlas
echo ""
echo "Step 3: Re-running L8 ATAC integration..."
python "${PROJECT_ROOT}/Analysis/ATAC/Integration/scripts/35_atac_integration.py"
echo "  L8 integration done"

# Step 4: Re-generate Fig 2 with corrected chromVAR
echo ""
echo "Step 4: Regenerating Fig 2..."
micromamba activate rnaseq
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
Rscript "${PROJECT_ROOT}/scripts/figures/fig2_compact.R" && echo "  Fig 2 SUCCEEDED" || echo "  Fig 2 FAILED"

echo ""
echo "============================================================"
echo "chromVAR post-processing COMPLETE: $(date)"
echo "============================================================"
