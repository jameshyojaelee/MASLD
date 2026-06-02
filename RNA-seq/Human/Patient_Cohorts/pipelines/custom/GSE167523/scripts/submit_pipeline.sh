#!/bin/bash
#SBATCH --job-name=liver_pipeline_GSE167523
#SBATCH --output=logs/submit_%j.out
#SBATCH --error=logs/submit_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00

set -euo pipefail

# USAGE: sbatch submit_pipeline.sh

# 1. Setup Environment
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

# Resolve Project Root (Robust to SLURM spooling)
# SLURM_SUBMIT_DIR points to the directory where sbatch was executed (scripts/)
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
PROJECT_ROOT=$(dirname "$SUBMIT_DIR")
cd "$PROJECT_ROOT"
echo "Running in PROJECT_ROOT: $PROJECT_ROOT"
# Infer dataset ID from path
DATASET_ID=$(basename $PROJECT_ROOT)
RESULTS_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/${DATASET_ID}"
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/${DATASET_ID}/fastq"
SAMPLESHEET="${PROJECT_ROOT}/metadata/samples.tsv"
SAMPLESHEET_WITH_PATHS="${PROJECT_ROOT}/metadata/samples_with_paths.tsv"

# 2. Check FASTQs and Update Samplesheet
if [[ -d "$FASTQ_DIR" && "$(ls -A $FASTQ_DIR)" ]]; then
    echo "Files found in $FASTQ_DIR. Updating samplesheet..."
    python scripts/generate_samplesheet.py "$SAMPLESHEET" "$FASTQ_DIR" "$SAMPLESHEET_WITH_PATHS"
else
    echo "Warning: No FASTQ files found. Pipeline may fail."
fi

# 3. Submit Snakemake (Master Process)
# We run snakemake here, which submits usage jobs to SLURM via --cluster
# Critical: Ensure --cluster argument uses 'cpu' partition if io is down

mkdir -p logs

echo "Starting Snakemake orchestrator..."
snakemake \
  -s workflow/Snakefile \
  --jobs 100 \
  --latency-wait 120 \
  --rerun-incomplete \
  --cluster-config workflow/cluster_config.yaml \
  --cluster "sbatch -J liver_{rule} -A {cluster.account} -p cpu --time {cluster.time} --mem {cluster.mem} {cluster.extra}" \
  --jobscript workflow/jobscript.sh \
  --printshellcmds

echo "Pipeline finished."
