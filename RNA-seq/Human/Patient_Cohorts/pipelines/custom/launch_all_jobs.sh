#!/bin/bash
set -euo pipefail

# Launch script for Human RNA-seq pipelines
# Submits PARALLEL ARRAY download jobs and then chains the pipeline submission

BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom"

DATASETS=("GSE126848" "GSE167523" "PRJNA512027" "GSE213621")
CHUNK_SIZE=5

for DATASET in "${DATASETS[@]}"; do
    echo "Processing $DATASET..."
    SCRIPT_DIR="${BASE_DIR}/${DATASET}/scripts"
    
    if [[ ! -d "$SCRIPT_DIR" ]]; then
        echo "Error: Directory $SCRIPT_DIR not found."
        continue
    fi
    
    cd "$SCRIPT_DIR"
    
    if [[ ! -f "download_fastq.sh" ]] || [[ ! -f "accession_list.txt" ]]; then
         echo "Error: Missing scripts or accession list in $SCRIPT_DIR"
         continue
    fi
    
    # Calculate Array Index Limit
    TOTAL_LINES=$(grep -cve '^\s*$' accession_list.txt || true) # Count non-empty lines
    if [[ "$TOTAL_LINES" -eq 0 ]]; then
        echo "Error: accession_list.txt is empty."
        continue
    fi
    
    # Calculate ceil(TOTAL_LINES / CHUNK_SIZE)
    ARRAY_LIMIT=$(( (TOTAL_LINES + CHUNK_SIZE - 1) / CHUNK_SIZE ))
    
    echo "  Found $TOTAL_LINES samples."
    echo "  Submitting Array Job 1-${ARRAY_LIMIT} (Chunk size $CHUNK_SIZE)..."
    
    # Submit Array Download Job
    # We explicitly pass the array config here
    DOWNLOAD_JOB_ID=$(sbatch --parsable --array=1-${ARRAY_LIMIT} download_fastq.sh)
    echo "  Submitted download job array: $DOWNLOAD_JOB_ID"
    
    # Create a wrapper sbatch script to launch the pipeline
    # We use quoted EOF to prevent variable expansion by the parent shell
    cat > launch_pipeline_wrapper.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=liver_launcher
#SBATCH --output=launcher_%j.out
#SBATCH --error=launcher_%j.err
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

echo "Pipeline launcher started."
# Ensure MICROMAMBA is available
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
if [ -x "$MICROMAMBA" ]; then
  eval "$($MICROMAMBA shell hook --shell bash)"
  micromamba activate rnaseq
else
  echo "Error: Micromamba not found at $MICROMAMBA"
  exit 1
fi

# Run the submit script
echo "Running submit_pipeline.sh..."
bash submit_pipeline.sh
EOF

    # Submit Pipeline Launcher with dependency
    # afterok:JOB_ID ensures it waits for ALL tasks in the array to succeed
    LAUNCHER_JOB_ID=$(sbatch --parsable --dependency=afterok:${DOWNLOAD_JOB_ID} launch_pipeline_wrapper.sh)
    echo "  Submitted pipeline launcher job: $LAUNCHER_JOB_ID (depends on $DOWNLOAD_JOB_ID)"
    echo "-----------------------------------"
done

echo "All jobs submitted."
