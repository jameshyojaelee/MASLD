#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=cpu,bigmem,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/quant_mouse_%A_%a.out
#SBATCH --error=logs/quant_mouse_%A_%a.err
#SBATCH --array=1-464

# Phase 1 — clean mouse transcript-level requant (vM38 primary index), per-dataset
# strandedness (strand_map.tsv from mouse_strand_infer). NB: -b 100 is a NO-OP on the
# kallisto/0.51.1 module (no HDF5) — point estimates only. Resume on abundance.tsv.
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
module load kallisto/0.51.1
IDX=$PROJ/data/reference/kallisto/gencode.vM38.kallisto.idx
DRIVER=mouse_driver.tsv
STRAND_MAP=strand_map_mouse.tsv     # cols: dataset  strand(unstranded|rf|fr)
OUT_ROOT=$PROJ/RNA-seq/results/isoform_diversity/mouse/kallisto

ROW=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+1' "$DRIVER")
[ -z "$ROW" ] && { echo "no driver row $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
DS=$(echo "$ROW" | cut -f1); SID=$(echo "$ROW" | cut -f2); LAY=$(echo "$ROW" | cut -f3)
R1=$(echo "$ROW" | cut -f4); R2=$(echo "$ROW" | cut -f5)
STR=$(awk -F'\t' -v d="$DS" '$1==d{print $2; exit}' "$STRAND_MAP"); STR=${STR:-unstranded}
FLAG=""; [ "$STR" = rf ] && FLAG="--rf-stranded"; [ "$STR" = fr ] && FLAG="--fr-stranded"

OUT_DIR="$OUT_ROOT/$DS/$SID"; mkdir -p "$OUT_DIR"
[ -s "$OUT_DIR/abundance.tsv" ] && { echo "$DS/$SID exists -- skip"; exit 0; }   # .h5 never written on this build
echo "[$(date)] $DS/$SID layout=$LAY strand=$STR"
if [ "$LAY" = "PAIRED" ] && [ -n "$R2" ]; then
  kallisto quant -i "$IDX" -o "$OUT_DIR" -t "$SLURM_CPUS_PER_TASK" -b 100 $FLAG "$R1" "$R2"
else
  kallisto quant -i "$IDX" -o "$OUT_DIR" -t "$SLURM_CPUS_PER_TASK" -b 100 $FLAG --single -l 200 -s 30 "$R1"
fi
python3 -c "import json;d=json.load(open('$OUT_DIR/run_info.json'));print('p_pseudoaligned',d['p_pseudoaligned'])"
echo "[$(date)] DONE $OUT_DIR"
