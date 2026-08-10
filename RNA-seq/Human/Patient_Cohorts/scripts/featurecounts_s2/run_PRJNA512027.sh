#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --partition=cpu,io,bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/featurecounts_s2/PRJNA512027_%j.out
#SBATCH --error=logs/featurecounts_s2/PRJNA512027_%j.err

echo "ERROR: PRJNA512027 is retired/noncanonical (L0/S0 library-preparation confound); recount is prohibited." >&2
exit 64

set -eo pipefail

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

set -u

DATASET="PRJNA512027"
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
BAM_DIR="${BASE}/results/${DATASET}/alignments/star"
OUT_DIR="${BASE}/results/${DATASET}/counts/featurecounts"
OUT_FILE="${OUT_DIR}/gene_counts.txt"

echo "=== featureCounts re-run: ${DATASET} (PE, -s 2) ==="
echo "Started: $(date)"

# Back up existing counts
if [[ -f "${OUT_FILE}" ]]; then
    cp "${OUT_FILE}" "${OUT_DIR}/gene_counts_s0_backup.txt"
    echo "Backed up existing gene_counts.txt -> gene_counts_s0_backup.txt"
fi
if [[ -f "${OUT_FILE}.summary" ]]; then
    cp "${OUT_FILE}.summary" "${OUT_DIR}/gene_counts_s0_backup.txt.summary"
fi

# Collect BAM files
BAMS=$(find "${BAM_DIR}" -name "*.Aligned.sortedByCoord.out.bam" | sort)
N_BAMS=$(echo "${BAMS}" | wc -l)
echo "Found ${N_BAMS} BAM files"

# Run featureCounts — PE mode (-p -B -s 2 reverse-stranded)
featureCounts \
    -T 8 \
    -p -B \
    -s 2 \
    -a "${GTF}" \
    -o "${OUT_FILE}" \
    ${BAMS}

# Verify output
if [[ -f "${OUT_FILE}" ]] && [[ $(stat -c%s "${OUT_FILE}") -gt 0 ]]; then
    echo "SUCCESS: gene_counts.txt written ($(stat -c%s "${OUT_FILE}") bytes)"
    echo "Gene count: $(tail -n +3 "${OUT_FILE}" | wc -l)"
else
    echo "ERROR: Output file missing or empty"
    exit 1
fi

echo "Finished: $(date)"
