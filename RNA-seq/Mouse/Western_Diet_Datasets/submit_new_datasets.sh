#!/bin/bash
# =============================================================================
# Chain: Kallisto quant -> tximport -> DE for GSE292565 + GSE246328
# =============================================================================
# Prerequisites: FASTQ downloads must be complete for both datasets.
# This script submits Kallisto quant as dependent on download jobs if provided.
#
# Usage:
#   # After downloads complete (or to check and submit):
#   bash submit_new_datasets.sh
#
#   # Chain after download jobs:
#   bash submit_new_datasets.sh 16555957 16555962
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WDDIR="${BASE}/RNA-seq/Mouse/Western_Diet_Datasets"

DL_JOB_292565="${1:-}"
DL_JOB_246328="${2:-}"

echo "============================================="
echo "  New WD Datasets Pipeline"
echo "============================================="
echo "  $(date)"
echo ""

# ---- Step 1: Kallisto quant ----
echo "--- Step 1: Kallisto Quantification ---"

# GSE292565 (paired-end, 24 samples)
N292=$(wc -l < "${WDDIR}/GSE292565/metadata/sample_list.txt")
DEP292=""
[[ -n "${DL_JOB_292565}" ]] && DEP292="--dependency=afterok:${DL_JOB_292565}"

KQ292=$(sbatch ${DEP292} \
    --array=1-${N292} \
    --job-name=kall_GSE292565 \
    --partition=cpu,io \
    --qos=nslab \
    "${WDDIR}/run_kallisto_quant.sh" GSE292565 \
    | grep -oP '\d+')
echo "  GSE292565: Kallisto job ${KQ292} (${N292} paired-end samples)"

# GSE246328 (single-end, 55 samples)
N246=$(wc -l < "${WDDIR}/GSE246328/metadata/sample_list.txt")
DEP246=""
[[ -n "${DL_JOB_246328}" ]] && DEP246="--dependency=afterok:${DL_JOB_246328}"

KQ246=$(sbatch ${DEP246} \
    --array=1-${N246} \
    --job-name=kall_GSE246328 \
    --partition=cpu,io \
    --qos=nslab \
    "${WDDIR}/run_kallisto_quant.sh" GSE246328 --single -l 200 -s 30 \
    | grep -oP '\d+')
echo "  GSE246328: Kallisto job ${KQ246} (${N246} single-end samples)"

# ---- Step 2: Per-dataset tximport + DE ----
echo ""
echo "--- Step 2: Per-Dataset tximport + DE ---"

DE292=$(sbatch \
    --dependency=afterok:${KQ292} \
    --job-name=de_GSE292565 \
    "${WDDIR}/submit_tximport_de.sh" GSE292565 \
    | grep -oP '\d+')
echo "  GSE292565: DE job ${DE292} (depends on Kallisto ${KQ292})"

DE246=$(sbatch \
    --dependency=afterok:${KQ246} \
    --job-name=de_GSE246328 \
    "${WDDIR}/submit_tximport_de.sh" GSE246328 \
    | grep -oP '\d+')
echo "  GSE246328: DE job ${DE246} (depends on Kallisto ${KQ246})"

# ---- Step 3: Combined tximport + unified DE ----
echo ""
echo "--- Step 3: Combined tximport + Unified DE ---"
echo "  After per-dataset DE completes, rebuild combined txi:"
echo "    sbatch --dependency=afterok:${DE292}:${DE246} \\"
echo "      --partition=cpu --mem=32G --time=48:00:00 \\"
echo "      --job-name=txi_rebuild \\"
echo "      --wrap='eval \"\$(micromamba shell hook -s bash)\"; micromamba activate rnaseq; Rscript ${WDDIR}/build_combined_txi.R && Rscript ${WDDIR}/run_wd_kallisto_de.R'"

echo ""
echo "============================================="
echo "  Job Summary"
echo "============================================="
echo "  GSE292565 download: ${DL_JOB_292565:-already_done}"
echo "  GSE246328 download: ${DL_JOB_246328:-already_done}"
echo "  GSE292565 Kallisto: ${KQ292}"
echo "  GSE246328 Kallisto: ${KQ246}"
echo "  GSE292565 DE:       ${DE292}"
echo "  GSE246328 DE:       ${DE246}"
echo ""
echo "Monitor: squeue -u \$(whoami) --name=kall_GSE292565,kall_GSE246328,de_GSE292565,de_GSE246328"
echo "============================================="
