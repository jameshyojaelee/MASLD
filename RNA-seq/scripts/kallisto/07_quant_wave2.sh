#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/quant_w2_%A_%a.out
#SBATCH --error=logs/quant_w2_%A_%a.err

# Wave 2 - Kallisto quant for 4 disease-only staging cohorts (423 samples).
# Reuses index from wave 1 (worktree). Outputs to main repo's kallisto dir.
# Driver: sample_driver_wave2.tsv (cohort, sample_id, layout, fastq_r1, fastq_r2)
# Submit: sbatch --array=1-423%4 07_quant_wave2.sh

set -eo pipefail
module load kallisto/0.51.1
set -u

IDX=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/results/kallisto/_index/gencode.v49.kidx
DRIVER=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto/sample_driver_wave2.tsv
OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/kallisto

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

# Skip if already done (abundance.tsv exists and is non-empty)
if [[ -s "$OUT_DIR/abundance.tsv" ]]; then
  echo "[$(date)] $COHORT/$SAMPLE already quantified -- skipping"
  exit 0
fi

# Verify FASTQ exists
if [[ ! -f "$R1" ]]; then
  echo "[ERROR] R1 FASTQ not found: $R1" >&2; exit 1
fi

echo "[$(date)] Cohort=$COHORT Sample=$SAMPLE Layout=$LAYOUT"
echo "  R1=$R1"
echo "  R2=$R2"
echo "  OUT=$OUT_DIR"

if [[ "$LAYOUT" == "PAIRED" && -n "$R2" && "$R2" != "NA" ]]; then
  if [[ ! -f "$R2" ]]; then
    echo "[ERROR] R2 FASTQ not found: $R2" >&2; exit 1
  fi
  kallisto quant -i "$IDX" -o "$OUT_DIR" -t "$SLURM_CPUS_PER_TASK" -b 30 "$R1" "$R2"
else
  # Single-end: fragment length 200 +/- 30 (standard Illumina short-insert)
  kallisto quant -i "$IDX" -o "$OUT_DIR" -t "$SLURM_CPUS_PER_TASK" -b 30 \
    --single -l 200 -s 30 "$R1"
fi

echo "[$(date)] Done: $COHORT/$SAMPLE"
