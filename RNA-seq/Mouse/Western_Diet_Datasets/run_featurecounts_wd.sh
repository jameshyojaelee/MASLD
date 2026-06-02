#!/bin/bash
#SBATCH --job-name=fCounts_WD
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/featurecounts_wd_%j.out
#SBATCH --error=logs/featurecounts_wd_%j.err
# Submit after STAR: sbatch --dependency=afterok:JOB1:JOB2:... run_featurecounts_wd.sh

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

PROJECT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WD_DIR="${PROJECT}/RNA-seq/Mouse/Western_Diet_Datasets"
GTF="${PROJECT}/RNA-seq/Mouse/Public_Diet_Models/reference/raw/gencode.vM38.annotation.gtf"
THREADS=${SLURM_CPUS_PER_TASK:-16}

for DATASET_ID in GSE220575 GSE246088 GSE246328 GSE292565 GSE305484; do
    STAR_DIR="${WD_DIR}/${DATASET_ID}/alignments/star"
    OUT_DIR="${WD_DIR}/${DATASET_ID}/counts/featurecounts"
    mkdir -p "${OUT_DIR}"

    BAMS=$(find "${STAR_DIR}" -name "*.Aligned.sortedByCoord.out.bam" | sort)
    NBAMS=$(echo "${BAMS}" | wc -l)
    echo "=== ${DATASET_ID}: ${NBAMS} BAMs ==="

    if [[ "${NBAMS}" -eq 0 ]]; then
        echo "SKIPPING ${DATASET_ID}: no BAMs found" >&2
        continue
    fi

    SAMPLE_LIST="${WD_DIR}/${DATASET_ID}/metadata/sample_list.txt"
    NCOLS=$(head -1 "${SAMPLE_LIST}" | awk -F'\t' '{print NF}')

    if [[ "${NCOLS}" -ge 3 ]]; then
        echo "PE mode (-p --countReadPairs)"
        featureCounts -T "${THREADS}" -s 0 -p --countReadPairs \
            -a "${GTF}" \
            -o "${OUT_DIR}/gene_counts.txt" \
            ${BAMS}
    else
        echo "SE mode (no -p flag)"
        featureCounts -T "${THREADS}" -s 0 \
            -a "${GTF}" \
            -o "${OUT_DIR}/gene_counts.txt" \
            ${BAMS}
    fi

    echo "${DATASET_ID} featureCounts DONE: $(wc -l < "${OUT_DIR}/gene_counts.txt") genes"
done

echo "=== All featureCounts complete ==="
