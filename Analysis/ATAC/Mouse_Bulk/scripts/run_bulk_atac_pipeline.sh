#!/bin/bash
#SBATCH --job-name=bulk_atac_pipeline
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Mouse_Bulk/logs/bulk_atac_pipeline_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Mouse_Bulk/logs/bulk_atac_pipeline_%j.err

# =============================================================================
# Module 1: Bulk Mouse ATAC-seq Pipeline (GSE246213)
# Runs the full Snakemake pipeline: fastp → Bowtie2 → filter → MACS2 →
# consensus peaks → featureCounts → QC
# =============================================================================

set -euo pipefail

# --- Paths (absolute to avoid SLURM working directory issues) ---
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_DIR="${PROJECT_ROOT}/Analysis/ATAC/Mouse_Bulk"

cd "$PIPELINE_DIR"

echo "============================================"
echo "Bulk Mouse ATAC-seq Pipeline"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "============================================"

# --- Load HPC modules ---
module purge
module load Bowtie2/2.5.4-linux-x86_64
module load SAMtools/1.21
module load BEDTools/2.31.0-GCC-12.3.0
module load MACS2/2.2.9.1-foss-2023b
module load picard/3.0.0-Java-17
module load deepTools/3.5.6
module load Subread/2.0.4-GCC-11.3.0
module load MultiQC/1.22.3-foss-2023b

# --- Ensure fastp is available ---
FASTP_BIN="$PIPELINE_DIR/bin/fastp"
if [[ ! -x "$FASTP_BIN" ]]; then
    echo "Installing fastp (standalone binary)..."
    mkdir -p "$PIPELINE_DIR/bin"
    wget -q -O "$FASTP_BIN" http://opengene.org/fastp/fastp
    chmod +x "$FASTP_BIN"
fi
export PATH="$PIPELINE_DIR/bin:$PATH"

# --- Verify tools ---
echo ""
echo "Tool versions:"
echo "  bowtie2:      $(bowtie2 --version 2>&1 | head -1)"
echo "  samtools:     $(samtools --version | head -1)"
echo "  bedtools:     $(bedtools --version)"
echo "  macs2:        $(macs2 --version 2>&1)"
echo "  picard:       $(picard MarkDuplicates --version 2>&1 | tail -1 || echo 'available')"
echo "  deeptools:    $(deeptools --version 2>&1)"
echo "  featureCounts: $(featureCounts -v 2>&1 | head -1)"
echo "  fastp:        $(fastp --version 2>&1)"
echo ""

# --- Verify reference genome ---
BT2_INDEX="/gpfs/commons/home/jameslee/reference_genome/bowtie2/GRCm39/GRCm39"
if [[ ! -f "${BT2_INDEX}.1.bt2" ]]; then
    echo "ERROR: Bowtie2 index not found at ${BT2_INDEX}"
    exit 1
fi
echo "Bowtie2 index: OK"

# --- Generate chrom sizes if missing ---
CHROM_SIZES="/gpfs/commons/home/jameslee/reference_genome/bowtie2/GRCm39/GRCm39.chrom.sizes"
if [[ ! -f "$CHROM_SIZES" ]]; then
    echo "Generating chrom.sizes from Bowtie2 index..."
    bowtie2-inspect -s "$BT2_INDEX" | grep "^Sequence" | awk '{print $2"\t"$3}' > "$CHROM_SIZES"
fi
echo "Chrom sizes: OK"

# --- Create output directories ---
RESULTS_DIR="$PIPELINE_DIR/results"
mkdir -p "$RESULTS_DIR"/{trimmed,alignment,peaks,qc}
mkdir -p "$PIPELINE_DIR/logs"

# --- Activate snakemake env ---
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo ""
echo "Snakemake version: $(snakemake --version)"
echo ""

# --- Dry run first ---
echo "=== DRY RUN ==="
snakemake --snakefile workflow/Snakefile \
    --configfile workflow/config.yaml \
    -n -p --reason 2>&1 | tail -20

echo ""
echo "=== STARTING PIPELINE ==="
echo ""

# --- Run pipeline ---
snakemake --snakefile workflow/Snakefile \
    --configfile workflow/config.yaml \
    --cores "$SLURM_CPUS_PER_TASK" \
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
echo "Consensus peaks: $(wc -l < "$RESULTS_DIR/consensus_peaks.bed" 2>/dev/null || echo 'not found') peaks"
echo "Count matrix: $(ls -lh "$RESULTS_DIR/consensus_peak_counts.txt" 2>/dev/null || echo 'not found')"
echo "QC report: $RESULTS_DIR/qc/multiqc_report.html"
echo ""
echo "Next step: Run DiffBind analysis"
echo "  Rscript scripts/03_diffbind_analysis.R"
