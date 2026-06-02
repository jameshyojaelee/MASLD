#!/bin/bash
# =============================================================================
# scaffold_new_human_dataset.sh
# Scaffold the directory structure for a new human RNA-seq dataset.
#
# Usage:
#   bash scripts/scaffold_new_human_dataset.sh GSE999999
#
# What this creates:
#   RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE999999/
#     metadata/           ← put SraRunTable.csv here
#     scripts/
#       accession_list.txt   ← one SRR ID per line
#       download_fastq.sh    ← SLURM array download script
#       submit_pipeline.sh   ← SLURM orchestrator
#       generate_samplesheet.py → symlink to shared script
#     workflow/
#       config.yaml          ← per-dataset config (edit project_root)
#       cluster_config.yaml  ← SLURM cluster config
#
# After running this script:
#   1. Copy SraRunTable.csv into metadata/
#   2. Fill accession_list.txt (one SRR per line)
#   3. Edit workflow/config.yaml (project_root path)
#   4. Add dataset entry to config/human_datasets.yaml
#   5. Add a harmonize() case to 00_harmonize_metadata.R
#   6. Submit: cd scripts && sbatch download_fastq.sh
# =============================================================================

set -euo pipefail

DATASET_ID="${1:-}"
if [[ -z "$DATASET_ID" ]]; then
    echo "Usage: $0 <DATASET_ID>   (e.g., GSE999999)"
    exit 1
fi

# Locate project root (script lives in scripts/ one level below root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
HUMAN_BASE="${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts"
SHARED_SCRIPTS="${HUMAN_BASE}/pipelines/custom/scripts"
DEST="${HUMAN_BASE}/pipelines/custom/${DATASET_ID}"

if [[ -d "$DEST" ]]; then
    echo "ERROR: Dataset directory already exists: $DEST"
    exit 1
fi

echo "Creating scaffold for: $DATASET_ID"
echo "Destination: $DEST"

# ---------------------------------------------------------------------------
# Create directory structure
# ---------------------------------------------------------------------------
mkdir -p "${DEST}/metadata"
mkdir -p "${DEST}/scripts"
mkdir -p "${DEST}/workflow"

# ---------------------------------------------------------------------------
# workflow/config.yaml
# ---------------------------------------------------------------------------
cat > "${DEST}/workflow/config.yaml" << YAML
# Per-dataset Snakemake config — adjust project_root below.
project_root: ${HUMAN_BASE}/results/${DATASET_ID}
samplesheet: metadata/samples_with_paths.tsv

reference:
  gtf: /gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz

indices:
  star: /gpfs/commons/home/jameslee/reference_genome/gencode_v49/star_index

qc:
  fastqc_dir: qc/fastqc
  multiqc_dir: qc/multiqc

alignment:
  star_out: alignments/star
  featurecounts_out: counts/featurecounts

threads:
  fastqc: 4
  star: 16
  featurecounts: 8

params:
  featurecounts:
    strandedness: 0   # 0=unstranded, 1=stranded, 2=reverse-stranded

pipeline_logs: logs
YAML

# ---------------------------------------------------------------------------
# workflow/cluster_config.yaml
# ---------------------------------------------------------------------------
cat > "${DEST}/workflow/cluster_config.yaml" << YAML
__default__:
  account: ""
  time: "2:00:00"
  mem: "32G"
  extra: ""

star:
  time: "4:00:00"
  mem: "64G"

featurecounts:
  time: "2:00:00"
  mem: "32G"
YAML

# ---------------------------------------------------------------------------
# scripts/accession_list.txt — placeholder
# ---------------------------------------------------------------------------
cat > "${DEST}/scripts/accession_list.txt" << 'TXT'
# Fill this file with one SRR accession per line.
# Example:
# SRR1234567
# SRR1234568
TXT

# ---------------------------------------------------------------------------
# scripts/download_fastq.sh — copy from existing dataset as template
# ---------------------------------------------------------------------------
TEMPLATE_DOWNLOAD="${HUMAN_BASE}/pipelines/custom/GSE126848/scripts/download_fastq.sh"
if [[ -f "$TEMPLATE_DOWNLOAD" ]]; then
    # Replace GSE126848 with new dataset ID in the copy
    sed "s/GSE126848/${DATASET_ID}/g" "$TEMPLATE_DOWNLOAD" \
        > "${DEST}/scripts/download_fastq.sh"
    chmod +x "${DEST}/scripts/download_fastq.sh"
    echo "  Created: scripts/download_fastq.sh (from template)"
else
    echo "  WARNING: template download_fastq.sh not found at $TEMPLATE_DOWNLOAD"
    touch "${DEST}/scripts/download_fastq.sh"
fi

# ---------------------------------------------------------------------------
# scripts/submit_pipeline.sh
# ---------------------------------------------------------------------------
cat > "${DEST}/scripts/submit_pipeline.sh" << SBATCH
#!/bin/bash
#SBATCH --job-name=pipeline_${DATASET_ID}
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "\$(\$MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

SUBMIT_DIR="\${SLURM_SUBMIT_DIR:-\$(pwd)}"
PROJECT_DIR="\$(dirname "\$SUBMIT_DIR")"

cd "\$PROJECT_DIR"

DATASET_ID="${DATASET_ID}"
FASTQ_DIR="${HUMAN_BASE}/data/raw/\${DATASET_ID}/fastq"
SAMPLESHEET="\${PROJECT_DIR}/metadata/samples.tsv"
SAMPLESHEET_PATHS="\${PROJECT_DIR}/metadata/samples_with_paths.tsv"

# Generate samplesheet with FASTQ paths
if [[ -d "\$FASTQ_DIR" && "\$(ls -A \$FASTQ_DIR)" ]]; then
    echo "Updating samplesheet with FASTQ paths..."
    python "${SHARED_SCRIPTS}/generate_samplesheet.py" \
        "\$SAMPLESHEET" "\$FASTQ_DIR" "\$SAMPLESHEET_PATHS"
else
    echo "WARNING: No FASTQ files in \$FASTQ_DIR — pipeline may fail."
fi

mkdir -p logs

echo "Starting Snakemake..."
snakemake \\
  -s "${HUMAN_BASE}/pipelines/shared/Snakefile" \\
  --configfile workflow/config.yaml \\
  --jobs 100 \\
  --latency-wait 120 \\
  --rerun-incomplete \\
  --cluster-config workflow/cluster_config.yaml \\
  --cluster "sbatch -p cpu --time {cluster.time} --mem {cluster.mem} {cluster.extra}" \\
  --printshellcmds

echo "Pipeline finished."
SBATCH
chmod +x "${DEST}/scripts/submit_pipeline.sh"

# ---------------------------------------------------------------------------
# Symlink to shared scripts
# ---------------------------------------------------------------------------
if [[ -f "${SHARED_SCRIPTS}/generate_samplesheet.py" ]]; then
    ln -sf "${SHARED_SCRIPTS}/generate_samplesheet.py" \
        "${DEST}/scripts/generate_samplesheet.py"
fi

# ---------------------------------------------------------------------------
# Print next steps
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo " Scaffold created: $DEST"
echo "============================================================"
echo ""
echo " NEXT STEPS:"
echo ""
echo "  1. Download SRA metadata table from GEO/SRA:"
echo "       cp SraRunTable.csv ${DEST}/metadata/"
echo ""
echo "  2. Fill accession_list.txt:"
echo "       awk -F',' 'NR>1{print \$1}' ${DEST}/metadata/SraRunTable.csv \\"
echo "           > ${DEST}/scripts/accession_list.txt"
echo ""
echo "  3. Edit workflow/config.yaml if needed (paths, strandedness)"
echo ""
echo "  4. Add dataset entry to config/human_datasets.yaml:"
echo "       Edit: ${PROJECT_ROOT}/config/human_datasets.yaml"
echo ""
echo "  5. Add harmonize() case to 00_harmonize_metadata.R"
echo ""
echo "  6. Submit download jobs:"
echo "       cd ${DEST}/scripts"
echo "       TOTAL=\$(wc -l < accession_list.txt)"
echo "       sbatch --array=1-\$(( (TOTAL + 4) / 5 )) download_fastq.sh"
echo ""
echo "  7. After download completes, submit pipeline:"
echo "       sbatch submit_pipeline.sh"
echo ""
echo "  8. After featureCounts finishes, re-run integration:"
echo "       cd ${HUMAN_BASE}/analysis/integration/scripts"
echo "       sbatch run_full_5cohort_slurm.sh"
echo "============================================================"
