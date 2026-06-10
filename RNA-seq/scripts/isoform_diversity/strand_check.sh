#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/strand_check_%j.out
#SBATCH --error=logs/strand_check_%j.err

# Confirm library strandedness before launching the 661-sample array.
# Reverse-stranded (dUTP) libraries: --rf-stranded keeps the reads, --fr-stranded
# rejects most -> p_pseudoaligned(rf) ~ unstranded >> p_pseudoaligned(fr).
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
module load kallisto/0.51.1

IDX=$PROJ/RNA-seq/results/isoform_diversity/_index/gencode.v49.primary.kidx
OUT=$PROJ/RNA-seq/results/isoform_diversity/_strandcheck
mkdir -p "$OUT"

# one sample per PE cohort
for line in \
  "GSE213621" "GSE135251" "GSE130970"; do
  COH=$line
  ROW=$(awk -F'\t' -v c="$COH" 'NR>1 && $1==c {print; exit}' human_driver.tsv)
  SID=$(echo "$ROW" | cut -f2); R1=$(echo "$ROW" | cut -f4); R2=$(echo "$ROW" | cut -f5)
  echo "=== $COH / $SID ==="
  for MODE in unstranded rf fr; do
    FLAG=""; [ "$MODE" = "rf" ] && FLAG="--rf-stranded"; [ "$MODE" = "fr" ] && FLAG="--fr-stranded"
    o="$OUT/${COH}_${SID}_${MODE}"; mkdir -p "$o"
    kallisto quant -i "$IDX" -o "$o" -t "$SLURM_CPUS_PER_TASK" $FLAG "$R1" "$R2" 2>/dev/null
    p=$(python3 -c "import json;print(json.load(open('$o/run_info.json'))['p_pseudoaligned'])")
    echo "  $MODE: p_pseudoaligned=$p"
  done
done
echo "[$(date)] strand check DONE"
