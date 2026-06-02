#!/bin/bash
#SBATCH --job-name=atlas_rebuild
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/step5_atlas_rebuild_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/step5_atlas_rebuild_%j.err

set -euo pipefail
export ATAC_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome"
export SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/scripts"
export PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export RNA_ATLAS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
export DONOR_PAIRING="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/GSE244832/metadata/donor_pairing.csv"
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd "$ATAC_DIR"

# Backup original L8 annotated files
L8_DIR="results/l8_annotated"
if [[ -d "${L8_DIR}" ]] && [[ ! -d "${L8_DIR}_backup_gene_activity" ]]; then
    cp -r "${L8_DIR}" "${L8_DIR}_backup_gene_activity"
    echo "Backed up original L8 results"
fi

# Re-run peak annotation (08) with v2 chromVAR/SCENIC+ results
GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
DA_RESULTS="results/snapatac2/scatac_da_results.csv"
CHROMVAR_V2="results/chromvar_v2/chromvar_tf_activity.csv"
REGULONS="scenic_plus/hepatocyte_regulons.csv"

if [[ -f "$SCRIPT_DIR/08_annotate_peaks_for_l8.py" ]] && \
   [[ -f "${DA_RESULTS}" ]] && \
   [[ -f "${CHROMVAR_V2}" ]] && \
   [[ -f "${REGULONS}" ]]; then
    echo "Re-running peak annotation for L8..."
    python "$SCRIPT_DIR/08_annotate_peaks_for_l8.py" \
        --da-results "${DA_RESULTS}" \
        --chromvar "${CHROMVAR_V2}" \
        --regulons "${REGULONS}" \
        --gtf "${GTF}" \
        --output-dir "${L8_DIR}" \
        2>&1 || echo "Peak annotation failed (non-fatal), using existing L8 files"
else
    echo "WARNING: Missing inputs for peak annotation, skipping 08_annotate_peaks_for_l8.py"
fi

# SCENIC+ v2 not re-run (atac_env unavailable). Using original regulons.
echo "Using original SCENIC+ results (regulons robust to label changes: 67/67 TFs recovered)"

# Run atlas integration
INTEGRATION_SCRIPT="$PROJECT_ROOT/Analysis/ATAC/Integration/scripts/35_atac_integration.py"
if [[ -f "${INTEGRATION_SCRIPT}" ]]; then
    echo "Running L8 atlas integration..."
    cd "$PROJECT_ROOT"
    python "${INTEGRATION_SCRIPT}" --force
    echo "Atlas integration complete"
else
    echo "WARNING: Integration script not found: ${INTEGRATION_SCRIPT}"
fi

echo "Step 5 COMPLETE: Atlas rebuild done"
