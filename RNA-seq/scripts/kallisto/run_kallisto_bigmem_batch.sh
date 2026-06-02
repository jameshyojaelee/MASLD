#!/bin/bash
#SBATCH --job-name=kallisto-bigmem
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/kallisto_bigmem_%j.out
#SBATCH --error=RNA-seq/logs/kallisto_bigmem_%j.err

# Run remaining kallisto quant tasks (7-423) on bigmem node in parallel batches.
# Skips any sample already quantified (idempotent via abundance.tsv check).
# 16 parallel processes at a time (4 CPUs each = 64 CPUs total).

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
module load kallisto/0.51.1

IDX=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/results/kallisto/_index/gencode.v49.kidx
DRIVER=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto/sample_driver_wave2.tsv
OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/kallisto

MAX_PARALLEL=32
COMPLETED=0
FAILED=0
SKIPPED=0

quant_one() {
    local task_id=$1
    local row=$(awk -v i="$task_id" 'NR==i+1' "$DRIVER")
    local cohort=$(echo "$row" | cut -f1)
    local sample=$(echo "$row" | cut -f2)
    local layout=$(echo "$row" | cut -f3)
    local r1=$(echo "$row" | cut -f4)
    local r2=$(echo "$row" | cut -f5)
    local outdir="$OUT_ROOT/$cohort/$sample"
    mkdir -p "$outdir"

    if [[ -s "$outdir/abundance.tsv" ]]; then
        echo "[SKIP] $cohort/$sample already done"
        return 0
    fi

    if [[ "$layout" == "SINGLE" ]]; then
        kallisto quant -i "$IDX" -o "$outdir" --single -l 200 -s 20 -t 2 "$r1" 2>"$outdir/kallisto.err"
    else
        kallisto quant -i "$IDX" -o "$outdir" -t 2 "$r1" "$r2" 2>"$outdir/kallisto.err"
    fi

    if [[ -s "$outdir/abundance.tsv" ]]; then
        echo "[DONE] $cohort/$sample"
    else
        echo "[FAIL] $cohort/$sample"
        return 1
    fi
}

echo "[$(date)] Starting kallisto bigmem batch: tasks 7-423 (417 samples), $MAX_PARALLEL parallel"

running=0
for task_id in $(seq 7 423); do
    quant_one "$task_id" &
    ((running++))
    if (( running >= MAX_PARALLEL )); then
        wait -n 2>/dev/null || true
        ((running--))
    fi
done
wait

echo "[$(date)] All samples processed"

# Count results
total=$(seq 7 423 | wc -l)
done_count=$(find "$OUT_ROOT" -name "abundance.tsv" -newer "$DRIVER" | wc -l)
echo "Completed: $done_count / $total new samples"
echo "Total abundance.tsv files: $(find "$OUT_ROOT" -name "abundance.tsv" | wc -l)"
