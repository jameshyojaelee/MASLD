#!/bin/bash
#SBATCH --job-name=liver_pipeline_GSE162694
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
PROJECT_DIR="$(dirname "$SUBMIT_DIR")"

cd "$PROJECT_DIR"

DATASET_ID="GSE162694"
FASTQ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/${DATASET_ID}/fastq"
SAMPLESHEET="${PROJECT_DIR}/metadata/samples.tsv"
SAMPLESHEET_PATHS="${PROJECT_DIR}/metadata/samples_with_paths.tsv"

# Generate samplesheet with FASTQ paths
if [[ -d "$FASTQ_DIR" && "$(ls -A $FASTQ_DIR)" ]]; then
    echo "Updating samplesheet with FASTQ paths..."
    python "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/scripts/generate_samplesheet.py"         "$SAMPLESHEET" "$FASTQ_DIR" "$SAMPLESHEET_PATHS"
else
    echo "WARNING: No FASTQ files in $FASTQ_DIR — pipeline may fail."
fi

mkdir -p logs

echo "Starting Snakemake..."
snakemake \
  -s "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/shared/Snakefile" \
  --configfile workflow/config.yaml \
  --jobs 100 \
  --latency-wait 120 \
  --rerun-incomplete \
  --cluster-config workflow/cluster_config.yaml \
  --cluster "sbatch -J liver_{rule} -p cpu --time {cluster.time} --mem {cluster.mem} {cluster.extra}" \
  --printshellcmds

echo "Pipeline finished."
