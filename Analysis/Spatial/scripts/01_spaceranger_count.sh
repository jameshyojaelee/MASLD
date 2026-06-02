#!/usr/bin/env bash
#SBATCH --job-name=spaceranger_%a
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=Analysis/Spatial/logs/spaceranger_%A_%a.out
#SBATCH --error=Analysis/Spatial/logs/spaceranger_%A_%a.err

##############################################################################
# 01_spaceranger_count.sh — SpaceRanger alignment for Visium spatial data
#
# Usage:
#   sbatch --array=0-14 01_spaceranger_count.sh GSE192741
#   sbatch --array=0-34 01_spaceranger_count.sh HRA007511
##############################################################################

set -euo pipefail

DATASET="${1:?Usage: sbatch --array=0-N 01_spaceranger_count.sh <GSE192741|HRA007511>}"

echo "=== SpaceRanger Count ==="
echo "Dataset: ${DATASET}"
echo "Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Date: $(date)"
echo ""

# Load SpaceRanger
module load SpaceRanger/2.0.0-GCC-11.2.0

# Project paths
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SPATIAL_ROOT="Analysis/Spatial"
REFDIR="${HOME}/reference_genome/cellranger/refdata-gex-GRCh38-2024-A"

# Read sample ID for this array index
SAMPLE_LIST="${SPATIAL_ROOT}/metadata/${DATASET}_samples.txt"
if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "ERROR: Sample list not found: ${SAMPLE_LIST}"
    echo "Run 00c_rename_fastqs.py first to generate sample lists."
    exit 1
fi

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_LIST}")
if [ -z "${SAMPLE}" ]; then
    echo "ERROR: No sample at index ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST}"
    exit 1
fi
echo "Sample: ${SAMPLE}"

# Input FASTQ directory (symlinks created by 00c)
FASTQ_DIR="${SPATIAL_ROOT}/spaceranger_input/${DATASET}/${SAMPLE}"
if [ ! -d "${FASTQ_DIR}" ]; then
    echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}"
    echo "Run 00c_rename_fastqs.py first."
    exit 1
fi

# Check if FASTQs exist
N_FASTQ=$(ls "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null | wc -l)
if [ "${N_FASTQ}" -eq 0 ]; then
    echo "ERROR: No FASTQ files in ${FASTQ_DIR}"
    exit 1
fi
echo "FASTQ files: ${N_FASTQ}"

# Output directory
OUTBASE="${SPATIAL_ROOT}/results/spaceranger/${DATASET}"
mkdir -p "${OUTBASE}"

# Skip if already complete
if [ -f "${OUTBASE}/${SAMPLE}/outs/filtered_feature_bc_matrix.h5" ]; then
    echo "SKIP: ${SAMPLE} already processed (filtered_feature_bc_matrix.h5 exists)"
    exit 0
fi

# Check for histology image
IMAGE_DIR="${SPATIAL_ROOT}/metadata/histology/${DATASET}"
IMAGE_FLAG=""
if [ -f "${IMAGE_DIR}/${SAMPLE}.tif" ]; then
    IMAGE_FLAG="--image=${IMAGE_DIR}/${SAMPLE}.tif"
    echo "H&E image: ${IMAGE_DIR}/${SAMPLE}.tif"
elif [ -f "${IMAGE_DIR}/${SAMPLE}.jpg" ]; then
    IMAGE_FLAG="--image=${IMAGE_DIR}/${SAMPLE}.jpg"
    echo "H&E image: ${IMAGE_DIR}/${SAMPLE}.jpg"
elif [ -f "${IMAGE_DIR}/${SAMPLE}.png" ]; then
    IMAGE_FLAG="--image=${IMAGE_DIR}/${SAMPLE}.png"
    echo "H&E image: ${IMAGE_DIR}/${SAMPLE}.png"
else
    echo "WARN: No H&E image for ${SAMPLE}, using --unknown-slide visium-1"
    IMAGE_FLAG="--unknown-slide visium-1"
fi

# Run SpaceRanger count
cd "${OUTBASE}"

spaceranger count \
    --id="${SAMPLE}" \
    --transcriptome="${REFDIR}" \
    --fastqs="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/${FASTQ_DIR}" \
    ${IMAGE_FLAG} \
    --localcores=${SLURM_CPUS_PER_TASK} \
    --localmem=120

echo ""
echo "=== SpaceRanger Complete ==="
echo "Sample: ${SAMPLE}"
echo "Output: ${OUTBASE}/${SAMPLE}/outs/"
echo "End: $(date)"
