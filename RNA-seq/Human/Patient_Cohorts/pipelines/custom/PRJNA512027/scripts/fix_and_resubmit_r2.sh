#!/bin/bash
#SBATCH --job-name=liver_fix_prjna_r2
#SBATCH --output=logs/fix_prjna_r2_%j.out
#SBATCH --error=logs/fix_prjna_r2_%j.err
#SBATCH --partition=cpu
#SBATCH --time=06:00:00
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
STAR_DIR="${PROJECT_ROOT}/results/PRJNA512027/alignments/star"
PIPELINE_DIR="${PROJECT_ROOT}/pipelines/custom/PRJNA512027"

cd "$PROJECT_ROOT"

# 7 samples with corrupted FASTQs
SAMPLES=(SRR8378431 SRR8378451 SRR8378481 SRR8378547 SRR8378557 SRR8378560 SRR8378611)

# ── Step 1: Re-download all corrupted FASTQs ─────────────────
echo "=== Step 1: Re-downloading ${#SAMPLES[@]} corrupted samples ==="
for s in "${SAMPLES[@]}"; do
    echo "--- Processing $s ---"
    rm -f "${FASTQ_DIR}/${s}"*.fastq.gz "${FASTQ_DIR}/${s}"*.fastq
    rm -rf "${STAR_DIR}/${s}"
    
    if ! fasterq-dump --outdir "${FASTQ_DIR}" --split-files --threads 8 --progress "$s"; then
        echo "ERROR: fasterq-dump failed for $s — skipping"
        continue
    fi
    gzip -f "${FASTQ_DIR}/${s}"*.fastq
    
    count=$(ls "${FASTQ_DIR}/${s}"_*.fastq.gz 2>/dev/null | wc -l)
    echo "$s: $count FASTQs"
    if [ "$count" -ne 2 ]; then
        echo "WARNING: Expected 2 FASTQs for $s, got $count"
    fi
done

# ── Step 2: Verify integrity ─────────────────────────────────
echo "=== Step 2: Verifying FASTQ integrity ==="
PASS=0; FAIL=0
for s in "${SAMPLES[@]}"; do
    for f in "${FASTQ_DIR}/${s}"_*.fastq.gz; do
        if [ ! -f "$f" ]; then continue; fi
        if gzip -t "$f" 2>/dev/null; then
            PASS=$((PASS+1))
        else
            echo "CORRUPT: $f"
            FAIL=$((FAIL+1))
        fi
    done
done
echo "Integrity check: $PASS passed, $FAIL failed"

# ── Step 3: Unlock Snakemake & re-submit ─────────────────────
echo "=== Step 3: Unlocking Snakemake ==="
micromamba activate rnaseq
cd "${PIPELINE_DIR}"
snakemake --unlock -s workflow/Snakefile --configfile workflow/config.yaml 2>&1 || true

echo "=== Step 4: Re-submitting pipeline ==="
cd scripts
sbatch submit_pipeline.sh

echo "=== Fix round 2 complete ==="
