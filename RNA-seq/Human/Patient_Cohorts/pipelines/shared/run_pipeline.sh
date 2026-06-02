#!/bin/bash
#SBATCH --job-name=liver_pipeline
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1

# =============================================================================
# Master per-dataset RNA-seq pipeline launcher
#
# Usage (from a dataset dir):
#   sbatch ../../shared/run_pipeline.sh
#
# Or with a specific dataset dir:
#   sbatch --export=DATASET_DIR=/path/to/GSE999999 ../../shared/run_pipeline.sh
#
# Expects:
#   workflow/config.yaml          — dataset-specific config
#   metadata/samples_with_paths.tsv — sample manifest with FASTQ paths
#
# Pipeline steps (all parallelized via SLURM array-style cluster submission):
#   1. FastQC (per-sample, parallel)
#   2. STAR alignment (per-sample, parallel, 16 CPU / 64G each)
#   3. featureCounts (per-sample, parallel, 4 CPU / 16G each)
#   4. Merge featureCounts (single job, after all per-sample jobs)
#   5. MultiQC (single job, after merge)
# =============================================================================

set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

# Determine dataset directory
DATASET_DIR="${DATASET_DIR:-$(pwd)}"
cd "${DATASET_DIR}"

DATASET_ID=$(basename "${DATASET_DIR}")
echo "=== Pipeline: ${DATASET_ID} ==="
echo "Directory: ${DATASET_DIR}"
echo "Started: $(date)"

# Update job name for squeue readability
scontrol update JobId="${SLURM_JOB_ID}" JobName="liver_pipeline_${DATASET_ID}" 2>/dev/null || true

mkdir -p logs/cluster

# Use shared configs; fall back to local if shared doesn't exist
SHARED_DIR="$(dirname "$(dirname "${DATASET_DIR}")")/shared"
SNAKEFILE="${SHARED_DIR}/Snakefile"
CLUSTER_CONFIG="${SHARED_DIR}/cluster_config.yaml"

# Allow local cluster config override
if [ -f "workflow/cluster_config.yaml" ]; then
    CLUSTER_CONFIG="workflow/cluster_config.yaml"
    echo "Using local cluster config: ${CLUSTER_CONFIG}"
fi

echo "Snakefile: ${SNAKEFILE}"
echo "Cluster config: ${CLUSTER_CONFIG}"

snakemake \
    -s "${SNAKEFILE}" \
    --configfile workflow/config.yaml \
    --cluster-config "${CLUSTER_CONFIG}" \
    --jobs 50 \
    --max-jobs-per-second 2 \
    --max-status-checks-per-second 10 \
    --cluster "sbatch \
                --partition={cluster.partition} \
                --cpus-per-task={cluster.cpusPerTask} \
                --mem={cluster.mem} \
                --time={cluster.time} \
                --job-name=liver_{cluster.jobName} \
                --output={cluster.logOut} \
                --error={cluster.logErr} \
                --parsable" \
    --latency-wait 120 \
    --keep-going \
    --rerun-incomplete

echo "=== Pipeline complete: ${DATASET_ID} ==="
echo "Finished: $(date)"
