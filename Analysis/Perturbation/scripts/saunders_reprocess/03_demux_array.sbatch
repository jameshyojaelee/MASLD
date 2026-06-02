#!/bin/bash
#SBATCH --job-name=prefetch
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --array=1-32%4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Perturbation/scripts/saunders_reprocess/logs/03_demux_%A_%a.log

# Array job:
#   For each of 32 SRRs in GSM8478327 Perturb-seq:
#     1) prefetch SRA archive
#     2) fasterq-dump to FASTQ
#     3) Run demux_saunders_sgrna.py → per-SRR counts.tsv.gz + stats.json
#     4) Delete FASTQ + .sra to save disk
#
# Uses --qos=interactive but cap concurrency at 4 (the interactive limit) via %4.

set -euo pipefail
module load SRA-Toolkit/3.2.1

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
WORK=$PROJ/Analysis/Perturbation/scripts/saunders_reprocess
DATA=$PROJ/data/perturbation/datasets/saunders2025
MANIFEST=$WORK/srr_manifest.tsv

# Read SRR for this array index
SRR=$(awk -v idx="$SLURM_ARRAY_TASK_ID" '$1==idx {print $2}' "$MANIFEST")
if [ -z "$SRR" ]; then echo "No SRR for index $SLURM_ARRAY_TASK_ID"; exit 1; fi

SCRATCH=/scratch/jameslee/saunders_demux/$SRR
mkdir -p "$SCRATCH" "$DATA/processed/per_srr_counts"

OUT_FILE=$DATA/processed/per_srr_counts/${SRR}.counts.tsv.gz
if [ -f "$OUT_FILE" ]; then
    echo "[$(date +%T)] $SRR output exists — skipping"
    exit 0
fi

echo "[$(date +%T)] === Task $SLURM_ARRAY_TASK_ID: $SRR ==="

# Step 1: prefetch (skip if already there)
SRA_FILE=$DATA/raw/sra/$SRR/$SRR.sra
if [ ! -f "$SRA_FILE" ]; then
    echo "[$(date +%T)] prefetch $SRR"
    prefetch --output-directory "$DATA/raw/sra" --max-size 50g "$SRR"
fi

# Step 2: fasterq-dump (to scratch, not data dir — these files are large)
R1=$SCRATCH/${SRR}_1.fastq
R2=$SCRATCH/${SRR}_2.fastq
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "[$(date +%T)] fasterq-dump $SRR"
    fasterq-dump --split-files --threads "${SLURM_CPUS_PER_TASK:-4}" \
        --temp "$SCRATCH/tmp" --outdir "$SCRATCH" "$SRA_FILE"
fi

# Step 3: demux
echo "[$(date +%T)] demux $SRR"
python3 "$WORK/demux_saunders_sgrna.py" \
    --r1 "$R1" --r2 "$R2" --srr "$SRR" \
    --out-dir "$DATA/processed/per_srr_counts"

# Step 4: cleanup scratch FASTQ (keep .sra so we can re-run demux without re-downloading)
echo "[$(date +%T)] cleanup $SCRATCH"
rm -f "$R1" "$R2"
rm -rf "$SCRATCH/tmp"

echo "[$(date +%T)] === DONE: $SRR ==="
