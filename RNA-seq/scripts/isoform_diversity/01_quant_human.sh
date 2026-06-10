#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=cpu,bigmem,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/quant_human_%A_%a.out
#SBATCH --error=logs/quant_human_%A_%a.err
#SBATCH --array=1-661

# Phase 1 — clean human transcript-level requant for isoform DTU.
# Reverse-stranded (dUTP) PE libraries. NB: -b 100 is a NO-OP on the kallisto/0.51.1
# module (built without HDF5) — no abundance.h5, no bootstraps; point estimates only
# (satuRn+stageR is the complete DTU method). Resume on abundance.tsv (the file that IS written).
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
module load kallisto/0.51.1

IDX=$PROJ/RNA-seq/results/isoform_diversity/_index/gencode.v49.primary.kidx
DRIVER=human_driver.tsv
OUT_ROOT=$PROJ/RNA-seq/results/isoform_diversity/human/kallisto

ROW=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+1' "$DRIVER")
[ -z "$ROW" ] && { echo "[ERROR] no driver row for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
COHORT=$(echo "$ROW" | cut -f1); SAMPLE=$(echo "$ROW" | cut -f2)
R1=$(echo "$ROW" | cut -f4);    R2=$(echo "$ROW" | cut -f5)

OUT_DIR="$OUT_ROOT/$COHORT/$SAMPLE"; mkdir -p "$OUT_DIR"
# resume on abundance.tsv: this kallisto build never writes abundance.h5, so the old
# .h5 guard never fired and every re-run silently re-quantified all 661 samples.
if [[ -s "$OUT_DIR/abundance.tsv" ]]; then
  echo "[$(date)] $COHORT/$SAMPLE already quantified -- skipping"; exit 0
fi

echo "[$(date)] $COHORT/$SAMPLE  R1=$R1"
kallisto quant -i "$IDX" -o "$OUT_DIR" -t "$SLURM_CPUS_PER_TASK" \
  -b 100 --rf-stranded "$R1" "$R2"
# sanity: pseudoalignment rate
python3 -c "import json;d=json.load(open('$OUT_DIR/run_info.json'));print('p_pseudoaligned',d['p_pseudoaligned'],'n_targets',d['n_targets'])"
echo "[$(date)] DONE $OUT_DIR"
