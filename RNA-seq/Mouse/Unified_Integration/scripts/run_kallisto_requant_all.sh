#!/bin/bash
# =============================================================================
# Orchestrator: Kallisto re-quantification of all 9 existing mouse datasets
# =============================================================================
# Submits SLURM array jobs for each dataset, using the generic Kallisto quant
# script. PE datasets get direct quant; SE datasets get --single -l 200 -s 30.
#
# All outputs go to <dataset_dir>/quant/kallisto/<sample_id>/abundance.{tsv,h5}
# Existing STAR counts are NOT touched.
#
# Usage:
#   bash run_kallisto_requant_all.sh          # submit all 9
#   bash run_kallisto_requant_all.sh --dry    # print commands without submitting
# =============================================================================

set -euo pipefail

DRY_RUN=0
[[ "${1:-}" == "--dry" ]] && DRY_RUN=1

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MOUSE="$BASE/RNA-seq/Mouse"
SCRIPT="$MOUSE/Unified_Integration/scripts/run_kallisto_quant_generic.sh"
LOGDIR="$MOUSE/Unified_Integration/logs"
mkdir -p "$LOGDIR"

# Track submitted job IDs for downstream dependency
declare -a JOB_IDS=()

submit_job() {
    local DATASET="$1"
    local SAMPLE_LIST="$2"
    local OUTDIR="$3"
    local LAYOUT="$4"        # PE or SE
    local N_SAMPLES
    N_SAMPLES=$(wc -l < "$SAMPLE_LIST")

    if [[ "$N_SAMPLES" -eq 0 ]]; then
        echo "[SKIP] $DATASET: empty sample list"
        return
    fi

    mkdir -p "$OUTDIR"

    local EXTRA_ARGS=""
    local SE_NOTE=""
    if [[ "$LAYOUT" == "SE" ]]; then
        EXTRA_ARGS="--single -l 200 -s 30"
        SE_NOTE=" (single-end)"
    fi

    local CMD="sbatch --job-name=kq_${DATASET} \
        --partition=io --qos=nslab \
        --mem=16G --cpus-per-task=8 --time=48:00:00 \
        --output=${LOGDIR}/kq_${DATASET}_%A_%a.out \
        --error=${LOGDIR}/kq_${DATASET}_%A_%a.err \
        --array=1-${N_SAMPLES} \
        ${SCRIPT} ${SAMPLE_LIST} ${OUTDIR} ${EXTRA_ARGS}"

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[DRY] $DATASET: $N_SAMPLES samples${SE_NOTE}"
        echo "      $CMD"
    else
        echo "[SUBMIT] $DATASET: $N_SAMPLES samples${SE_NOTE}"
        JOB_ID=$(eval "$CMD" | grep -o '[0-9]*')
        JOB_IDS+=("$JOB_ID")
        echo "         Job ID: $JOB_ID (array 1-$N_SAMPLES)"
    fi
}

echo "============================================================"
echo "  Kallisto Re-quantification — 9 Existing Mouse Datasets"
echo "  $(date)"
echo "============================================================"
echo ""

# -----------------------------------------------------------------------
# 1. InHouse_MCD (PE, 12 samples)
# -----------------------------------------------------------------------
submit_job "InHouse_MCD" \
    "$MOUSE/InHouse_MCD/metadata/sample_list_kallisto.txt" \
    "$MOUSE/InHouse_MCD/quant/kallisto" \
    "PE"

# -----------------------------------------------------------------------
# 2. GSE156918 — Public_MCD (SE, 11 samples)
# -----------------------------------------------------------------------
submit_job "GSE156918" \
    "$MOUSE/Public_MCD/GSE156918/metadata/sample_list_kallisto.txt" \
    "$MOUSE/Public_MCD/GSE156918/quant/kallisto" \
    "SE"

# -----------------------------------------------------------------------
# 3. GSE205974 — Public_MCD (PE, 6 samples)
# -----------------------------------------------------------------------
submit_job "GSE205974" \
    "$MOUSE/Public_MCD/GSE205974/metadata/sample_list_kallisto.txt" \
    "$MOUSE/Public_MCD/GSE205974/quant/kallisto" \
    "PE"

# -----------------------------------------------------------------------
# 4. GSE159911 — LIDPAD (PE, 152 samples)
# -----------------------------------------------------------------------
submit_job "GSE159911" \
    "$MOUSE/Public_Diet_Models/GSE159911/metadata/sample_list_kallisto.txt" \
    "$MOUSE/Public_Diet_Models/GSE159911/quant/kallisto" \
    "PE"

# -----------------------------------------------------------------------
# 5. GSE162876 — CDAHFD + FPC (PE, 213 samples)
# -----------------------------------------------------------------------
submit_job "GSE162876" \
    "$MOUSE/Public_Diet_Models/GSE162876/metadata/sample_list_kallisto.txt" \
    "$MOUSE/Public_Diet_Models/GSE162876/quant/kallisto" \
    "PE"

# -----------------------------------------------------------------------
# 6. GSE224069 — HFD (SE, 43 samples)
# -----------------------------------------------------------------------
submit_job "GSE224069" \
    "$MOUSE/Public_Diet_Models/GSE224069/metadata/sample_list_kallisto.txt" \
    "$MOUSE/Public_Diet_Models/GSE224069/quant/kallisto" \
    "SE"

# -----------------------------------------------------------------------
# 7. GSE225616 — GAN/DIO-NASH (PE, 19 samples)
# -----------------------------------------------------------------------
submit_job "GSE225616" \
    "$MOUSE/Public_Diet_Models/GSE225616/metadata/sample_list_kallisto.txt" \
    "$MOUSE/Public_Diet_Models/GSE225616/quant/kallisto" \
    "PE"

# -----------------------------------------------------------------------
# 8. GSE263273 — AMLN_ob (SE, 28 samples)
# -----------------------------------------------------------------------
submit_job "GSE263273" \
    "$MOUSE/Public_Diet_Models/GSE263273/metadata/sample_list_kallisto.txt" \
    "$MOUSE/Public_Diet_Models/GSE263273/quant/kallisto" \
    "SE"

# -----------------------------------------------------------------------
# 9. GSE274914 — HFD timepoints (PE, 25 samples)
# -----------------------------------------------------------------------
submit_job "GSE274914" \
    "$MOUSE/Public_Diet_Models/GSE274914/metadata/sample_list_kallisto.txt" \
    "$MOUSE/Public_Diet_Models/GSE274914/quant/kallisto" \
    "PE"

echo ""
echo "============================================================"
if [[ $DRY_RUN -eq 1 ]]; then
    echo "  DRY RUN complete. Re-run without --dry to submit."
else
    echo "  Submitted ${#JOB_IDS[@]} array jobs."
    echo "  Job IDs: ${JOB_IDS[*]}"
    echo ""
    echo "  Monitor: squeue -u \$USER --name='kq_*'"
    echo "  Total samples: ~509 across 9 datasets"
fi
echo "============================================================"
