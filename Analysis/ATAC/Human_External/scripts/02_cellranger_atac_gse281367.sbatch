#!/bin/bash
# 02_cellranger_atac_gse281367.sbatch  --  Layer A A2 (alignment arm)
# Run cellranger-atac count on each GSE281367 donor (true 10x snATAC).
# Launch AFTER 01_download completes (all 12 donors have R1/R2/R3 fastqs).
#
#SBATCH --job-name=cellranger
#SBATCH --partition=cpu
#SBATCH --array=1-12
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/cr_%A_%a.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/cr_%A_%a.err

set -eo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
REF=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-atac-GRCh38-2024-A
module load cellranger-atac/2.1.0

META="$PROJ/data/GSE281367/metadata/donor_pairing.csv"
# DONOR may be passed via --export (per-donor dependency jobs); else derive from array id
if [ -z "$DONOR" ]; then
  ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$META")
  DONOR=$(echo "$ROW" | cut -d, -f1)
fi
FASTQS="$PROJ/data/raw/GSE281367/fastq/$DONOR"
WORK="$PROJ/Analysis/ATAC/Human_External/cellranger"
mkdir -p "$WORK"; cd "$WORK"

echo "[$(date)] cellranger-atac count donor=$DONOR host=$(hostname)"
[ -f "$FASTQS/${DONOR}_S1_L001_R2_001.fastq.gz" ] || { echo "ERROR: fastqs missing for $DONOR" >&2; exit 2; }
if [ -f "$WORK/$DONOR/outs/singlecell.csv" ]; then echo "$DONOR already done; skipping."; exit 0; fi
rm -rf "$WORK/$DONOR" 2>/dev/null || true   # clear any partial run

cellranger-atac count \
  --id="$DONOR" \
  --reference="$REF" \
  --fastqs="$FASTQS" \
  --sample="$DONOR" \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=120

echo "[$(date)] $DONOR done. Summary:"
grep -m1 . "$WORK/$DONOR/outs/summary.csv" 2>/dev/null || true
