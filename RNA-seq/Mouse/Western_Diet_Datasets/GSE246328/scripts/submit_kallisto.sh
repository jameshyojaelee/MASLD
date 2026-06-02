#!/bin/bash
# Submit Kallisto quantification for GSE246328 after downloads complete
# Usage: bash submit_kallisto.sh [download_job_id]
#   If download_job_id is provided, jobs will depend on it completing successfully.
#   NOTE: This dataset is SINGLE-END (NextSeq 500, 75 cycles).
#   Kallisto needs --single -l <frag_len> -s <frag_sd>.
#   Using l=200 s=30 as reasonable defaults for NEBNext Ultra II RNA-seq.

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets"
DATASET="GSE246328"
N_SAMPLES=$(wc -l < "${BASE}/${DATASET}/metadata/sample_list.txt")

echo "Submitting Kallisto quant for ${DATASET}: ${N_SAMPLES} single-end samples"

DEPEND_FLAG=""
if [[ ${1:-} ]]; then
    DEPEND_FLAG="--dependency=afterok:${1}"
    echo "  Dependency: afterok:${1}"
fi

sbatch ${DEPEND_FLAG} \
    --array=1-${N_SAMPLES} \
    --job-name=kall_${DATASET} \
    --partition=cpu,io \
    --qos=nslab \
    ${BASE}/run_kallisto_quant.sh ${DATASET} --single -l 200 -s 30

echo "Submitted."
