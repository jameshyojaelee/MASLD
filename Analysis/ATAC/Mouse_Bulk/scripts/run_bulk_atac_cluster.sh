#!/bin/bash
#SBATCH --job-name=atac_controller
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Mouse_Bulk/logs/atac_controller_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Mouse_Bulk/logs/atac_controller_%j.err

# =============================================================================
# Module 1: Bulk Mouse ATAC-seq Pipeline — CLUSTER MODE
# Lightweight controller that submits per-rule SLURM jobs via Snakemake.
# All 12 samples process in parallel on separate nodes (~12x speedup).
# =============================================================================

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_DIR="${PROJECT_ROOT}/Analysis/ATAC/Mouse_Bulk"

cd "$PIPELINE_DIR"

echo "============================================"
echo "Bulk Mouse ATAC-seq Pipeline (CLUSTER MODE)"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "============================================"

# --- Ensure fastp is available ---
FASTP_BIN="$PIPELINE_DIR/bin/fastp"
if [[ ! -x "$FASTP_BIN" ]]; then
    echo "Installing fastp (standalone binary)..."
    mkdir -p "$PIPELINE_DIR/bin"
    wget -q -O "$FASTP_BIN" http://opengene.org/fastp/fastp
    chmod +x "$FASTP_BIN"
fi

# --- Create log directories ---
mkdir -p logs/slurm

# --- Activate snakemake env ---
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Snakemake version: $(snakemake --version)"
echo ""

# --- Clean up incomplete outputs from previous run ---
echo "Cleaning incomplete outputs from cancelled job..."
# Remove any .bam.tmp files from partial bowtie2 runs
find results/alignment -name "*.tmp*" -delete 2>/dev/null || true
# Unlock stale locks from previously cancelled jobs
snakemake --snakefile workflow/Snakefile --configfile workflow/config.yaml --unlock 2>/dev/null || true
echo ""

# --- Run pipeline in cluster mode ---
snakemake --snakefile workflow/Snakefile \
    --configfile workflow/config.yaml \
    --cluster "sbatch \
        --partition=cpu \
        --cpus-per-task={threads} \
        --mem={resources.mem_mb}M \
        --time={resources.time} \
        --job-name=atac_{rule} \
        --output=$PIPELINE_DIR/logs/slurm/%j_{rule}_{wildcards}.out \
        --error=$PIPELINE_DIR/logs/slurm/%j_{rule}_{wildcards}.err" \
    --jobscript workflow/slurm_jobscript.sh \
    --jobs 12 \
    --latency-wait 120 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    2>&1

echo ""
echo "============================================"
echo "Pipeline completed: $(date)"
echo "============================================"

# --- Summary ---
echo ""
echo "=== OUTPUT SUMMARY ==="
echo "Consensus peaks: $(wc -l < "$PIPELINE_DIR/results/consensus_peaks.bed" 2>/dev/null || echo 'not found') peaks"
echo "Count matrix: $(ls -lh "$PIPELINE_DIR/results/consensus_peak_counts.txt" 2>/dev/null || echo 'not found')"
echo "QC report: $PIPELINE_DIR/results/qc/multiqc_report.html"
echo ""
echo "Next step: Run DiffBind analysis"
echo "  Rscript scripts/03_diffbind_analysis.R"
