#!/bin/bash
#SBATCH --job-name=V2ref
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# V2 target screen: one-time reference preparation.
#   1. fetch the HISAT2 prebuilt GRCh38 genome index (untrusted download)
#   2. build a HISAT2 index of the unmasked target references
#   3. tile the target references and align the tiles to the human genome
#   4. apply the low-complexity and human-homology masks
#   5. build the HISAT2 index actually used by the screen, on the masked FASTA
#
# Usage: sbatch prepare_reference.sbatch <RESULT_DIR>
# <RESULT_DIR> must already contain the prepared target-reference FASTA.

set -euo pipefail

RESULT_DIR="${1:?result dir required}"
BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V2_viral_screen

# hisat2 and hisat2-build are wrapper scripts that call `python` and `perl`
# from PATH, so the environment's bin directory has to come first.
export PATH="$BIN:$PATH"

DL="$RESULT_DIR/downloads_untrusted"
REFS="$RESULT_DIR/refs"
mkdir -p "$DL/human_hisat2_index" "$REFS"

VIRAL_FA="$DL/viral_refs/viral_refs.fasta"
test -s "$VIRAL_FA"

echo "[$(date -u +%FT%TZ)] step 1: human HISAT2 index"
HUMAN_PREFIX="$DL/human_hisat2_index/grch38/genome"
if [ ! -s "${HUMAN_PREFIX}.1.ht2" ]; then
  curl -sS -o "$DL/human_hisat2_index/grch38_genome.tar.gz" \
    "https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"
  sha256sum "$DL/human_hisat2_index/grch38_genome.tar.gz" \
    > "$REFS/human_hisat2_index.sha256"
  tar -xzf "$DL/human_hisat2_index/grch38_genome.tar.gz" \
    -C "$DL/human_hisat2_index"
fi
ls -la "$(dirname "$HUMAN_PREFIX")"

echo "[$(date -u +%FT%TZ)] step 2: unmasked target index"
"$BIN/hisat2-build" -p 4 "$VIRAL_FA" "$REFS/viral_unmasked" > "$REFS/hisat2_build_unmasked.log" 2>&1

echo "[$(date -u +%FT%TZ)] step 3: tile and align tiles to human"
"$BIN/python" -I "$SCRIPTS/mask_viral_reference.py" tiles \
  "$VIRAL_FA" "$REFS/viral_tiles.fasta"
"$BIN/hisat2" -p 4 -f -x "$HUMAN_PREFIX" -U "$REFS/viral_tiles.fasta" \
  --no-spliced-alignment -k 1 \
  -S "$REFS/viral_tiles_vs_human.sam" \
  --summary-file "$REFS/viral_tiles_vs_human.summary.txt"

echo "[$(date -u +%FT%TZ)] step 4: apply masks"
"$BIN/python" -I "$SCRIPTS/mask_viral_reference.py" apply \
  "$VIRAL_FA" "$REFS/viral_tiles_vs_human.sam" \
  "$REFS/viral_refs.masked.fasta" "$REFS/mask_report.tsv"
gzip -f "$REFS/viral_tiles_vs_human.sam"

echo "[$(date -u +%FT%TZ)] step 5: masked target index"
"$BIN/hisat2-build" -p 4 "$REFS/viral_refs.masked.fasta" \
  "$REFS/viral_masked" > "$REFS/hisat2_build_masked.log" 2>&1

"$BIN/samtools" faidx "$REFS/viral_refs.masked.fasta"
cat "$REFS/mask_report.tsv"
echo "[$(date -u +%FT%TZ)] done"
