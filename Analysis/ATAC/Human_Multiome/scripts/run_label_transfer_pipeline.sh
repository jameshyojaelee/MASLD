#!/bin/bash
# =============================================================================
# SLURM Orchestrator: scATAC Label Transfer Pipeline
#
# Chains 6 steps to replace gene-activity-based cell type annotation with
# donor-matched snRNA label transfer, then re-runs all downstream analyses.
#
# Steps:
#   1. Label transfer (02b) — project ATAC gene activity → RNA PCA → KNN transfer
#   2. Peak re-calling (02c) — MACS3 with corrected cell types
#   3. chromVAR re-run (03) — TF motif enrichment with corrected labels
#   4. SCENIC+ re-run (04) — hepatocyte GRN with corrected hepatocyte selection
#   5. Atlas rebuild (35) — update S5 epigenomic columns
#   6. Sensitivity analysis (05) — compare old vs new downstream results
#
# Usage:
#   cd Analysis/ATAC/Human_Multiome
#   bash scripts/run_label_transfer_pipeline.sh
#
# Total wall time: ~20h (chained via SLURM dependencies)
# =============================================================================

set -euo pipefail

# --- Configuration ---
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATAC_DIR="${PROJECT_ROOT}/Analysis/ATAC/Human_Multiome"
SCRIPT_DIR="${ATAC_DIR}/scripts"
LOG_DIR="${ATAC_DIR}/logs"
LABEL_TRANSFER_DIR="${ATAC_DIR}/results/label_transfer"
SBATCH_DIR="${ATAC_DIR}/scripts/sbatch_label_transfer"

RNA_ATLAS="${PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
DONOR_PAIRING="${PROJECT_ROOT}/data/GSE244832/metadata/donor_pairing.csv"

mkdir -p "${LOG_DIR}" "${LABEL_TRANSFER_DIR}" "${SBATCH_DIR}"

echo "============================================================"
echo "scATAC Label Transfer Pipeline"
echo "  Project:    ${PROJECT_ROOT}"
echo "  ATAC dir:   ${ATAC_DIR}"
echo "  RNA atlas:  ${RNA_ATLAS}"
echo "  Started:    $(date)"
echo "============================================================"

# --- Validate inputs ---
FAIL=0
for F in "${RNA_ATLAS}" "${DONOR_PAIRING}"; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: Required input not found: $F"
        FAIL=1
    fi
done
if [[ ! -d "${ATAC_DIR}/results/snapatac2/per_donor" ]]; then
    echo "ERROR: Per-donor h5ad directory not found: ${ATAC_DIR}/results/snapatac2/per_donor"
    FAIL=1
fi
if [[ ! -f "${ATAC_DIR}/results/snapatac2/snapatac2_processed.h5ad" ]] && \
   [[ ! -f "${ATAC_DIR}/results/snapatac2/snapatac2_processed_fixed.h5ad" ]]; then
    echo "ERROR: Processed ATAC h5ad not found in ${ATAC_DIR}/results/snapatac2/"
    FAIL=1
fi
if [[ $FAIL -eq 1 ]]; then
    echo "Aborting: missing required inputs."
    exit 1
fi
echo "Input validation passed."

# ============================================================================
# Step 0: Extract RNA reference (requires rapids_singlecell env for anndata compat)
# ============================================================================
REF_OUTPUT="${LABEL_TRANSFER_DIR}/rna_reference_GSE244832.h5ad"
if [[ -f "${REF_OUTPUT}" ]]; then
    echo "Pre-extracted RNA reference found: ${REF_OUTPUT}"
    echo "  Skipping Step 0."
else
    echo "Extracting RNA reference from atlas (Step 0)..."
    cat > "${SBATCH_DIR}/step0_extract_ref.sbatch" << 'STEP0_EOF'
#!/bin/bash
#SBATCH --job-name=extract_ref
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=01:00:00

set -euo pipefail
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd "$ATAC_DIR"
python "$SCRIPT_DIR/02a_extract_rna_reference.py" \
    --atlas "$RNA_ATLAS" \
    --donor_pairing "$DONOR_PAIRING" \
    --output results/label_transfer/rna_reference_GSE244832.h5ad

echo "Step 0 COMPLETE: RNA reference extracted"
STEP0_EOF

    # Inject env vars and log paths
    sed -i "/^set -euo pipefail$/a\\
export ATAC_DIR=\"${ATAC_DIR}\"\\
export SCRIPT_DIR=\"${SCRIPT_DIR}\"\\
export PROJECT_ROOT=\"${PROJECT_ROOT}\"\\
export RNA_ATLAS=\"${RNA_ATLAS}\"\\
export DONOR_PAIRING=\"${DONOR_PAIRING}\"" "${SBATCH_DIR}/step0_extract_ref.sbatch"
    sed -i "s|^#SBATCH --time=.*|&\n#SBATCH --output=${LOG_DIR}/step0_extract_ref_%j.out\n#SBATCH --error=${LOG_DIR}/step0_extract_ref_%j.err|" "${SBATCH_DIR}/step0_extract_ref.sbatch"

    JOB0=$(sbatch --parsable "${SBATCH_DIR}/step0_extract_ref.sbatch")
    echo "Step 0 (Extract Reference): SLURM ${JOB0}"
    STEP1_DEP="--dependency=afterok:${JOB0}"
fi

# ============================================================================
# Generate per-step sbatch scripts (proper #!/bin/bash shebang)
# ============================================================================

# --- Step 1: Label Transfer ---
cat > "${SBATCH_DIR}/step1_label_transfer.sbatch" << 'STEP1_EOF'
#!/bin/bash
#SBATCH --job-name=label_transfer
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=08:00:00

set -euo pipefail
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd "$ATAC_DIR"
python "$SCRIPT_DIR/02b_label_transfer_from_rna.py" \
    --atac_dir results/snapatac2 \
    --rna_atlas "$RNA_ATLAS" \
    --donor_pairing "$DONOR_PAIRING" \
    --output_dir results/label_transfer \
    --k_neighbors 10 \
    --n_hvgs 3000 \
    --n_pcs 30

echo "Step 1 COMPLETE: Label transfer done"
STEP1_EOF

# --- Step 2: Peak Re-Calling ---
cat > "${SBATCH_DIR}/step2_peak_recall.sbatch" << 'STEP2_EOF'
#!/bin/bash
#SBATCH --job-name=peak_recall
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00

set -euo pipefail
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd "$ATAC_DIR"
python "$SCRIPT_DIR/02c_rerun_peaks_with_transferred_labels.py" \
    --label_transferred results/label_transfer/snapatac2_label_transferred.h5ad \
    --per_donor_dir results/snapatac2/per_donor \
    --old_peak_dir results/snapatac2/cell_type_peak_sets \
    --output_dir results/label_transfer

echo "Step 2 COMPLETE: Peak re-calling done"
STEP2_EOF

# --- Step 3: chromVAR v2 ---
cat > "${SBATCH_DIR}/step3_chromvar_v2.sbatch" << 'STEP3_EOF'
#!/bin/bash
#SBATCH --job-name=chromvar_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00

set -euo pipefail
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd "$ATAC_DIR"
RELABELED="results/label_transfer/snapatac2_relabeled.h5ad"
if [[ ! -f "${RELABELED}" ]]; then
    echo "ERROR: relabeled h5ad not found. Step 2 may have failed."
    exit 1
fi

python "$SCRIPT_DIR/03_chromvar_motifs.py" \
    --input "${RELABELED}" \
    --output-dir results/chromvar_v2 \
    --genome hg38

echo "Step 3 COMPLETE: chromVAR re-run done"
STEP3_EOF

# --- Step 4: SCENIC+ v2 ---
cat > "${SBATCH_DIR}/step4_scenic_v2.sbatch" << 'STEP4_EOF'
#!/bin/bash
#SBATCH --job-name=scenic_v2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00

set -euo pipefail
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate atac_env
CONDA_PREFIX_PATH="$(micromamba info | grep 'env location' | awk '{print $NF}')"
if [[ -d "${CONDA_PREFIX_PATH}/lib" ]]; then
    export LD_LIBRARY_PATH="${CONDA_PREFIX_PATH}/lib:${LD_LIBRARY_PATH:-}"
fi

cd "$ATAC_DIR"
RELABELED="results/label_transfer/snapatac2_relabeled.h5ad"
if [[ ! -f "${RELABELED}" ]]; then
    echo "ERROR: relabeled h5ad not found. Step 2 may have failed."
    exit 1
fi

python "$SCRIPT_DIR/04_scenic_plus_grn.py" \
    --atac-input "${RELABELED}" \
    --rna-dir cellranger_arc \
    --output-dir results/scenic_plus_v2 \
    --n-cells 15000 \
    --n-topics "10,20,30,40" \
    --n-iter 500 \
    --n-hvg 3000 \
    --genome hg38

echo "Step 4 COMPLETE: SCENIC+ re-run done"
STEP4_EOF

# --- Step 5: Atlas Rebuild ---
cat > "${SBATCH_DIR}/step5_atlas_rebuild.sbatch" << 'STEP5_EOF'
#!/bin/bash
#SBATCH --job-name=atlas_rebuild
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail
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
REGULONS_V2="results/scenic_plus_v2/hepatocyte_regulons.csv"

if [[ -f "$SCRIPT_DIR/08_annotate_peaks_for_l8.py" ]] && \
   [[ -f "${DA_RESULTS}" ]] && \
   [[ -f "${CHROMVAR_V2}" ]] && \
   [[ -f "${REGULONS_V2}" ]]; then
    echo "Re-running peak annotation for L8..."
    python "$SCRIPT_DIR/08_annotate_peaks_for_l8.py" \
        --da-results "${DA_RESULTS}" \
        --chromvar "${CHROMVAR_V2}" \
        --regulons "${REGULONS_V2}" \
        --gtf "${GTF}" \
        --output-dir "${L8_DIR}" \
        2>&1 || echo "Peak annotation failed (non-fatal), using existing L8 files"
else
    echo "WARNING: Missing inputs for peak annotation, skipping 08_annotate_peaks_for_l8.py"
fi

# Backup original SCENIC+ results
SCENIC_DIR="results/scenic_plus"
if [[ -d "${SCENIC_DIR}" ]] && [[ ! -d "${SCENIC_DIR}_backup_gene_activity" ]]; then
    cp -r "${SCENIC_DIR}" "${SCENIC_DIR}_backup_gene_activity"
fi

# Copy v2 SCENIC+ results over original location for atlas integration
if [[ -d "results/scenic_plus_v2" ]]; then
    for f in results/scenic_plus_v2/*.csv; do
        [ -f "$f" ] && cp "$f" "results/scenic_plus/$(basename "$f")"
    done
    echo "Copied SCENIC+ v2 results to scenic_plus/"
fi

# Run atlas integration
INTEGRATION_SCRIPT="$PROJECT_ROOT/Analysis/ATAC/Integration/scripts/35_atac_integration.py"
if [[ -f "${INTEGRATION_SCRIPT}" ]]; then
    echo "Running L8 atlas integration..."
    cd "$PROJECT_ROOT"
    python "${INTEGRATION_SCRIPT}"
    echo "Atlas integration complete"
else
    echo "WARNING: Integration script not found: ${INTEGRATION_SCRIPT}"
fi

echo "Step 5 COMPLETE: Atlas rebuild done"
STEP5_EOF

# --- Step 6: Sensitivity ---
cat > "${SBATCH_DIR}/step6_sensitivity.sbatch" << 'STEP6_EOF'
#!/bin/bash
#SBATCH --job-name=sensitivity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00

set -euo pipefail
module purge 2>/dev/null
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

cd "$ATAC_DIR"
python "$SCRIPT_DIR/05_annotation_sensitivity.py" \
    --label_transferred results/label_transfer/snapatac2_label_transferred.h5ad \
    --chromvar_old results/chromvar/chromvar_tf_activity.csv \
    --chromvar_new results/chromvar_v2/chromvar_tf_activity.csv \
    --scenic_old results/scenic_plus_backup_gene_activity/hepatocyte_regulons.csv \
    --scenic_new results/scenic_plus_v2/hepatocyte_regulons.csv \
    --output_dir results/label_transfer/sensitivity

echo "Step 6 COMPLETE: Sensitivity analysis done"
STEP6_EOF

# ============================================================================
# Inject environment variables into sbatch scripts and submit with dependencies
# ============================================================================

# All scripts need these variables — inject them after the #SBATCH block
for SCRIPT in "${SBATCH_DIR}"/step*.sbatch; do
    sed -i "/^set -euo pipefail$/a\\
export ATAC_DIR=\"${ATAC_DIR}\"\\
export SCRIPT_DIR=\"${SCRIPT_DIR}\"\\
export PROJECT_ROOT=\"${PROJECT_ROOT}\"\\
export RNA_ATLAS=\"${RNA_ATLAS}\"\\
export DONOR_PAIRING=\"${DONOR_PAIRING}\"" "$SCRIPT"
done

# Add log paths to each script
for SCRIPT in "${SBATCH_DIR}"/step*.sbatch; do
    STEP_NAME=$(basename "$SCRIPT" .sbatch)
    sed -i "s|^#SBATCH --time=.*|&\n#SBATCH --output=${LOG_DIR}/${STEP_NAME}_%j.out\n#SBATCH --error=${LOG_DIR}/${STEP_NAME}_%j.err|" "$SCRIPT"
done

# --- Submit with dependency chain ---
JOB1=$(sbatch --parsable ${STEP1_DEP:-} "${SBATCH_DIR}/step1_label_transfer.sbatch")
echo "Step 1 (Label Transfer): SLURM ${JOB1}"

JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} "${SBATCH_DIR}/step2_peak_recall.sbatch")
echo "Step 2 (Peak Re-Calling):  SLURM ${JOB2} (after ${JOB1})"

JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} "${SBATCH_DIR}/step3_chromvar_v2.sbatch")
echo "Step 3 (chromVAR v2):      SLURM ${JOB3} (after ${JOB2})"

JOB4=$(sbatch --parsable --dependency=afterok:${JOB2} "${SBATCH_DIR}/step4_scenic_v2.sbatch")
echo "Step 4 (SCENIC+ v2):       SLURM ${JOB4} (after ${JOB2})"

JOB5=$(sbatch --parsable --dependency=afterok:${JOB3}:${JOB4} "${SBATCH_DIR}/step5_atlas_rebuild.sbatch")
echo "Step 5 (Atlas Rebuild):    SLURM ${JOB5} (after ${JOB3},${JOB4})"

JOB6=$(sbatch --parsable --dependency=afterok:${JOB5} "${SBATCH_DIR}/step6_sensitivity.sbatch")
echo "Step 6 (Sensitivity):      SLURM ${JOB6} (after ${JOB5})"

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "============================================================"
echo "All 6 jobs submitted with SLURM dependencies:"
echo "  Step 1: ${JOB1}  Label Transfer"
echo "  Step 2: ${JOB2}  Peak Re-Calling       (after ${JOB1})"
echo "  Step 3: ${JOB3}  chromVAR v2           (after ${JOB2})"
echo "  Step 4: ${JOB4}  SCENIC+ v2            (after ${JOB2})"
echo "  Step 5: ${JOB5}  Atlas Rebuild          (after ${JOB3},${JOB4})"
echo "  Step 6: ${JOB6}  Sensitivity Analysis   (after ${JOB5})"
echo ""
echo "Dependency chain:"
echo "  1 -> 2 -> 3 \\"
echo "              5 -> 6"
echo "  1 -> 2 -> 4 /"
echo ""
echo "Monitor: squeue -u \$(whoami) --format='%.10i %.20j %.8T %.10M %.10l'"
echo "Logs:    ${LOG_DIR}/"
echo "Sbatch:  ${SBATCH_DIR}/"
echo "============================================================"
