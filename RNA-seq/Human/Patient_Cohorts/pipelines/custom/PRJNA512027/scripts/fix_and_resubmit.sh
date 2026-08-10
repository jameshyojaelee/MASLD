#!/bin/bash
#SBATCH --job-name=liver_fix_prjna
#SBATCH --output=logs/fix_prjna_%j.out
#SBATCH --error=logs/fix_prjna_%j.err
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=8
#SBATCH --account=nslab
#SBATCH --qos=nslab

echo "ERROR: PRJNA512027 is retired/noncanonical; resubmission is prohibited." >&2
exit 64

set -euo pipefail

# ── Environment ──────────────────────────────────────────────
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
eval "$($MICROMAMBA shell hook --shell bash)"
set +u
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl
set -u

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
FASTQ_DIR="${PROJECT_ROOT}/data/raw/PRJNA512027/fastq"
PIPELINE_DIR="${PROJECT_ROOT}/pipelines/custom/PRJNA512027"

cd "$PROJECT_ROOT"

# ── Step 1: Fix corrupted/missing FASTQs ─────────────────────
echo "=== Step 1: Re-downloading corrupted/missing FASTQs ==="

# SRR8378572: _2.fastq.gz is truncated
echo "Removing corrupted SRR8378572 files..."
rm -f "${FASTQ_DIR}/SRR8378572"*.fastq.gz "${FASTQ_DIR}/SRR8378572"*.fastq
echo "Re-downloading SRR8378572..."
fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --progress SRR8378572
gzip -f "${FASTQ_DIR}/SRR8378572"*.fastq
echo "SRR8378572 done."

# SRR8378555: _2.fastq.gz is missing
echo "Removing partial SRR8378555 files..."
rm -f "${FASTQ_DIR}/SRR8378555"*.fastq.gz "${FASTQ_DIR}/SRR8378555"*.fastq
echo "Re-downloading SRR8378555..."
fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --progress SRR8378555
gzip -f "${FASTQ_DIR}/SRR8378555"*.fastq
echo "SRR8378555 done."

# Verify
echo "=== Verifying FASTQs ==="
for s in SRR8378572 SRR8378555; do
    count=$(ls "${FASTQ_DIR}/${s}"_*.fastq.gz 2>/dev/null | wc -l)
    echo "$s: $count FASTQs"
    if [ "$count" -ne 2 ]; then
        echo "ERROR: Expected 2 FASTQs for $s, got $count"
        exit 1
    fi
done

# ── Step 2: Clean up failed STAR output ──────────────────────
echo "=== Step 2: Cleaning up failed STAR output for SRR8378572 ==="
rm -rf "${PROJECT_ROOT}/results/PRJNA512027/alignments/star/SRR8378572"

# ── Step 3: Unlock Snakemake ─────────────────────────────────
echo "=== Step 3: Unlocking Snakemake ==="
micromamba activate rnaseq
cd "${PIPELINE_DIR}"
snakemake --unlock -s workflow/Snakefile --configfile workflow/config.yaml 2>&1 || true

# ── Step 4: Re-submit pipeline ───────────────────────────────
echo "=== Step 4: Re-submitting pipeline ==="
cd scripts
sbatch submit_pipeline.sh

echo "=== Fix and resubmit complete ==="
