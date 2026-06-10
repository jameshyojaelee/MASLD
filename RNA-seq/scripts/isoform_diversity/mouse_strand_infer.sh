#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=cpu,bigmem,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/mouse_strand_%A_%a.out
#SBATCH --error=logs/mouse_strand_%A_%a.err
#SBATCH --array=1-9

# Infer per-dataset strandedness for the 9 mouse datasets (mixed library prep).
# One representative sample per dataset x {unstranded, rf, fr}; highest
# p_pseudoaligned wins. Uses the existing primary-clean vM38 index.
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
module load kallisto/0.51.1
IDX=$PROJ/data/reference/kallisto/gencode.vM38.kallisto.idx
OUT=$PROJ/RNA-seq/results/isoform_diversity/mouse/_strandcheck
mkdir -p "$OUT"

DS=$(awk -F'\t' 'NR>1{print $1}' mouse_driver.tsv | sort -u | awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i')
[ -z "$DS" ] && { echo "no dataset for task $SLURM_ARRAY_TASK_ID"; exit 0; }
ROW=$(awk -F'\t' -v d="$DS" 'NR>1 && $1==d{print; exit}' mouse_driver.tsv)
SID=$(echo "$ROW" | cut -f2); LAY=$(echo "$ROW" | cut -f3)
R1=$(echo "$ROW" | cut -f4); R2=$(echo "$ROW" | cut -f5)
echo "=== dataset=$DS sample=$SID layout=$LAY ==="
for MODE in unstranded rf fr; do
  FLAG=""; [ "$MODE" = rf ] && FLAG="--rf-stranded"; [ "$MODE" = fr ] && FLAG="--fr-stranded"
  o="$OUT/${DS}_${MODE}"; mkdir -p "$o"
  if [ "$LAY" = "PAIRED" ] && [ -n "$R2" ]; then
    kallisto quant -i "$IDX" -o "$o" -t "$SLURM_CPUS_PER_TASK" $FLAG "$R1" "$R2" 2>/dev/null
  else
    kallisto quant -i "$IDX" -o "$o" -t "$SLURM_CPUS_PER_TASK" $FLAG --single -l 200 -s 30 "$R1" 2>/dev/null
  fi
  p=$(python3 -c "import json;print(json.load(open('$o/run_info.json'))['p_pseudoaligned'])")
  echo "  $DS $MODE p_pseudoaligned=$p"
done
