#!/bin/bash
# =============================================================================
# Submit all Western Diet dataset pipelines: SRA download -> Kallisto quant
# =============================================================================
# This orchestrator submits SLURM jobs for each step:
#   1. SRA download (fasterq-dump) — all 3 datasets in parallel
#   2. Generate sample_list.txt (after downloads complete)
#   3. Kallisto quantification (after sample lists generated)
#
# Prerequisites:
#   - Kallisto index: data/reference/kallisto/gencode.vM38.kallisto.idx
#     (submit build_kallisto_index.sh first if missing)
#   - srr_list.txt in each dataset's metadata/ dir
#
# Usage:
#   cd RNA-seq/Mouse/Western_Diet_Datasets
#   bash submit_all_pipelines.sh
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WDDIR="${BASE}/RNA-seq/Mouse/Western_Diet_Datasets"
INDEX="${BASE}/data/reference/kallisto/gencode.vM38.kallisto.idx"

# Ensure log directory exists
mkdir -p "${WDDIR}/logs"

echo "============================================="
echo "  Western Diet Pipeline Orchestrator"
echo "============================================="
echo "  Date:      $(date)"
echo "  Base:      ${WDDIR}"
echo ""

# ---- Check Kallisto index ----
if [[ ! -s "${INDEX}" ]]; then
    echo "WARNING: Kallisto index missing or empty."
    echo "  Submit: sbatch data/reference/kallisto/build_kallisto_index.sh"
    echo "  Then re-run this script after the index is built."
    echo ""
    echo "  Proceeding with downloads only (Kallisto quant will be submitted manually)."
    INDEX_READY=0
else
    echo "Kallisto index: OK ($(du -h ${INDEX} | cut -f1))"
    INDEX_READY=1
fi
echo ""

# ---- Step 1: Submit SRA downloads ----
echo "--- Step 1: SRA Downloads ---"

declare -A DL_JOBS

for DS in GSE220575 GSE246088 GSE305484; do
    SRR_LIST="${WDDIR}/${DS}/metadata/srr_list.txt"
    N_SAMPLES=$(wc -l < "${SRR_LIST}")

    # Check how many are already downloaded
    FASTQDIR="${WDDIR}/${DS}/fastq"
    N_DONE=0
    while IFS= read -r SRR; do
        if [[ -f "${FASTQDIR}/${SRR}_1.fastq.gz" && -f "${FASTQDIR}/${SRR}_2.fastq.gz" ]]; then
            R1_SIZE=$(stat -c%s "${FASTQDIR}/${SRR}_1.fastq.gz" 2>/dev/null || echo 0)
            if [[ ${R1_SIZE} -gt 1000000 ]]; then
                N_DONE=$((N_DONE + 1))
            fi
        fi
    done < "${SRR_LIST}"

    if [[ ${N_DONE} -eq ${N_SAMPLES} ]]; then
        echo "  ${DS}: All ${N_SAMPLES} FASTQs already downloaded. Skipping."
        DL_JOBS[${DS}]="DONE"
    else
        REMAINING=$((N_SAMPLES - N_DONE))
        echo "  ${DS}: ${N_DONE}/${N_SAMPLES} done, submitting ${REMAINING} remaining..."

        JOB_ID=$(sbatch \
            --array=1-${N_SAMPLES} \
            --job-name="sra_${DS}" \
            --output="${WDDIR}/logs/sra_${DS}_%A_%a.out" \
            --error="${WDDIR}/logs/sra_${DS}_%A_%a.err" \
            "${WDDIR}/download_sra.sh" "${DS}" \
            | grep -oP '\d+')

        DL_JOBS[${DS}]="${JOB_ID}"
        echo "    Job ID: ${JOB_ID} (array 1-${N_SAMPLES})"
    fi
done

echo ""

# ---- Step 2: Submit Kallisto quant (dependent on downloads) ----
if [[ ${INDEX_READY} -eq 1 ]]; then
    echo "--- Step 2: Kallisto Quant (chained after downloads) ---"

    for DS in GSE220575 GSE246088 GSE305484; do
        N_SAMPLES=$(wc -l < "${WDDIR}/${DS}/metadata/srr_list.txt")

        if [[ "${DL_JOBS[${DS}]}" == "DONE" ]]; then
            # Downloads already done; generate sample list and submit quant directly
            echo "  ${DS}: Generating sample_list.txt..."
            bash "${WDDIR}/generate_sample_lists.sh" "${DS}"

            SAMPLE_LIST="${WDDIR}/${DS}/metadata/sample_list.txt"
            N_READY=$(wc -l < "${SAMPLE_LIST}")

            if [[ ${N_READY} -gt 0 ]]; then
                KQ_JOB=$(sbatch \
                    --array=1-${N_READY} \
                    --job-name="kallisto_${DS}" \
                    --output="${WDDIR}/logs/kallisto_${DS}_%A_%a.out" \
                    --error="${WDDIR}/logs/kallisto_${DS}_%A_%a.err" \
                    "${WDDIR}/run_kallisto_quant.sh" "${DS}" \
                    | grep -oP '\d+')
                echo "    Kallisto Job: ${KQ_JOB} (array 1-${N_READY})"
            fi
        else
            # Chain Kallisto after download completes
            echo "  ${DS}: Will chain after download job ${DL_JOBS[${DS}]}"
            echo "    NOTE: Run after downloads complete:"
            echo "      bash generate_sample_lists.sh ${DS}"
            echo "      N=\$(wc -l < ${DS}/metadata/sample_list.txt)"
            echo "      sbatch --array=1-\${N} run_kallisto_quant.sh ${DS}"
        fi
    done
else
    echo "--- Step 2: Skipped (Kallisto index not ready) ---"
    echo "  After index and downloads are complete, run:"
    echo "    bash generate_sample_lists.sh all"
    echo "    for DS in GSE220575 GSE246088 GSE305484; do"
    echo "      N=\$(wc -l < \${DS}/metadata/sample_list.txt)"
    echo "      sbatch --array=1-\${N} run_kallisto_quant.sh \${DS}"
    echo "    done"
fi

echo ""
echo "============================================="
echo "  Summary"
echo "============================================="
for DS in GSE220575 GSE246088 GSE305484; do
    echo "  ${DS}: download=${DL_JOBS[${DS}]}"
done
echo ""
echo "Monitor with:"
echo "  squeue -u \$(whoami) --name=sra_GSE220575,sra_GSE246088,sra_GSE305484,kallisto_GSE220575,kallisto_GSE246088,kallisto_GSE305484"
echo "============================================="
