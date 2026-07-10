#!/bin/bash
# run_nt_design.sh — build the mouse vM38 transcriptome Bowtie index (once) and
# design 500 non-targeting Cas13 control guides screened against it.
#   srun --partition=cpu --qos=interactive --ntasks=1 --cpus-per-task=8 --mem=64G \
#        --time=2:00:00 bash Cas13_Library_Design/scripts/guides/run_nt_design.sh
set -eo pipefail

ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
REF=/gpfs/commons/home/jameslee/reference_genome/refdata-cellranger-GRCm39-vM38
GENOME="$REF/fasta/genome.fa"
GTFGZ="$REF/genes/genes.gtf.gz"
CACHE="$ROOT/Cas13_Library_Design/data/guides/cache/nt_screen"
OUT="$ROOT/Cas13_Library_Design/data/guides/cas13_nt_controls_vM38.csv"
THREADS=8
mkdir -p "$CACHE"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq          # gffread + python/pandas
module load Bowtie/1.3.1-GCC-11.3.0  # bowtie + bowtie-build on PATH

# --- 1. transcriptome FASTA (all GENCODE vM38 transcripts) -------------------
if [ ! -s "$CACHE/transcriptome.fa" ]; then
  echo "=== [1/3] building transcriptome FASTA (gffread) ==="
  zcat "$GTFGZ" > "$CACHE/genes.gtf"
  gffread -w "$CACHE/transcriptome.fa" -g "$GENOME" "$CACHE/genes.gtf"
  rm -f "$CACHE/genes.gtf"
  echo "  transcripts: $(grep -c '^>' "$CACHE/transcriptome.fa")"
else echo "=== [1/3] transcriptome.fa cached ==="; fi

# --- 2. Bowtie index ---------------------------------------------------------
if [ ! -s "$CACHE/tx_idx.1.ebwt" ]; then
  echo "=== [2/3] bowtie-build transcriptome index ==="
  bowtie-build --threads "$THREADS" "$CACHE/transcriptome.fa" "$CACHE/tx_idx" >/dev/null
else echo "=== [2/3] bowtie index cached ==="; fi

# --- 3. generate + screen + emit 500 NT guides -------------------------------
echo "=== [3/3] design + transcriptome screen ==="
python "$ROOT/Cas13_Library_Design/scripts/guides/design_nt_controls.py" \
  --idx "$CACHE/tx_idx" --n 500 --out "$OUT" --workdir "$CACHE" --threads "$THREADS"
echo "=== DONE -> $OUT ==="
