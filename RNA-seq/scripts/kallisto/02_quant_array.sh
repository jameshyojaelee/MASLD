#!/bin/bash
#SBATCH --job-name=B1_kallisto_quant
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/quant_%A_%a.out
#SBATCH --error=logs/quant_%A_%a.err

# B1 - Kallisto quant array over 5 mega cohorts (860 samples pre-QC).
# Driver list: a unified TSV with columns cohort, sample_id, layout, fastq_r1, fastq_r2
# Strand: all 5 cohorts use featureCounts -s 0 (unstranded); kallisto default (no --rf/--fr) matches.
# Bootstrap=30 for downstream sleuth/tximport variance estimates.

set -euo pipefail
module load kallisto/0.51.1

WT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
IDX="$WT/RNA-seq/results/kallisto/_index/gencode.v49.kidx"
DRIVER="$WT/RNA-seq/scripts/kallisto/sample_driver.tsv"
OUT_ROOT="$WT/RNA-seq/results/kallisto"

# Pick row from driver (skip header). Array task ID = row number (1-based).
ROW=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+1' "$DRIVER")
if [[ -z "$ROW" ]]; then
  echo "[ERROR] No driver row for task $SLURM_ARRAY_TASK_ID" >&2; exit 1
fi
COHORT=$(echo "$ROW" | cut -f1)
SAMPLE=$(echo "$ROW" | cut -f2)
LAYOUT=$(echo "$ROW" | cut -f3)
R1=$(echo "$ROW" | cut -f4)
R2=$(echo "$ROW" | cut -f5)

OUT_DIR="$OUT_ROOT/$COHORT/$SAMPLE"
mkdir -p "$OUT_DIR"

if [[ -s "$OUT_DIR/abundance.h5" ]]; then
  echo "[$(date)] $COHORT/$SAMPLE already quantified -- skipping"
  exit 0
fi

echo "[$(date)] Cohort=$COHORT Sample=$SAMPLE Layout=$LAYOUT"
echo "  R1=$R1"
echo "  R2=$R2"

if [[ "$LAYOUT" == "PAIRED" && -n "$R2" && "$R2" != "NA" ]]; then
  kallisto quant -i "$IDX" -o "$OUT_DIR" -t "$SLURM_CPUS_PER_TASK" -b 30 "$R1" "$R2"
else
  # Single-end: kallisto requires --fragment-length and --sd. Use library-typical defaults.
  # GSE126848 / GSE162694 / PRJNA512027 are SE on HiSeq -> ~200 +/- 30 typical.
  kallisto quant -i "$IDX" -o "$OUT_DIR" -t "$SLURM_CPUS_PER_TASK" -b 30 \
    --single -l 200 -s 30 "$R1"
fi

echo "[$(date)] Done: $OUT_DIR"
