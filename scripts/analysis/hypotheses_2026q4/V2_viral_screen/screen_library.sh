#!/bin/bash
# V2 target screen: screen ONE library for target class B, target class C and target class E transcripts.
#
# Usage:
#   screen_library.sh <SRR> <COHORT> <RESULT_DIR> [FASTQ1] [FASTQ2]
#
# The STAR BAMs for these cohorts were produced without --outSAMunmapped, so
# they hold no unmapped reads. The screen therefore reads the FASTQ.
#
# Steps
#   1. align all reads end-to-end to the masked target references (HISAT2,
#      --no-spliced-alignment), keeping only aligned records
#   2. collapse the hits to unique read sequences
#   3. re-align those unique sequences to the human genome (HISAT2, spliced
#      alignment allowed) and drop any sequence the human genome explains at
#      least as well
#   4. apply the prespecified positivity rule
#
# Everything is written under <RESULT_DIR>/libraries/<COHORT>/<SRR>/.

set -euo pipefail

SRR="${1:?SRR required}"
COHORT="${2:?cohort required}"
RESULT_DIR="${3:?result dir required}"
FQ1_IN="${4:-}"
FQ2_IN="${5:-}"

BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V2_viral_screen
RAW=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# hisat2 is a wrapper script that calls `python` and `perl` from PATH.
export PATH="$BIN:$PATH"

REFS="$RESULT_DIR/refs"
VIRAL_IDX="$REFS/viral_masked"
HUMAN_IDX="$RESULT_DIR/downloads_untrusted/human_hisat2_index/grch38/genome"
ACC_TABLE="$RESULT_DIR/downloads_untrusted/viral_refs/viral_refs_accessions.tsv"
MASK_REPORT="$REFS/mask_report.tsv"

OUT="$RESULT_DIR/libraries/$COHORT/$SRR"
mkdir -p "$OUT"

# ---------------------------------------------------------------- input FASTQ
if [ -n "$FQ1_IN" ]; then
  FQ1="$FQ1_IN"; FQ2="$FQ2_IN"
else
  FQ1=""; FQ2=""
  if [ -s "$RAW/$COHORT/fastq/${SRR}_1.fastq.gz" ]; then
    FQ1="$RAW/$COHORT/fastq/${SRR}_1.fastq.gz"
    [ -s "$RAW/$COHORT/fastq/${SRR}_2.fastq.gz" ] && FQ2="$RAW/$COHORT/fastq/${SRR}_2.fastq.gz"
  elif [ -s "$RAW/$COHORT/fastq/${SRR}.fastq.gz" ]; then
    FQ1="$RAW/$COHORT/fastq/${SRR}.fastq.gz"
  fi
fi
test -s "$FQ1" || { echo "no FASTQ for $SRR in $COHORT" >&2; exit 3; }

if [ -n "$FQ2" ]; then
  LAYOUT=paired
  READ_ARGS=(-1 "$FQ1" -2 "$FQ2")
else
  LAYOUT=single
  READ_ARGS=(-U "$FQ1")
fi

{
  echo "library=$SRR"
  echo "cohort=$COHORT"
  echo "layout=$LAYOUT"
  echo "fastq1=$FQ1"
  echo "fastq2=${FQ2:-}"
  echo "threads=$THREADS"
  echo "started_utc=$(date -u +%FT%TZ)"
  echo "hisat2_version=$("$BIN/hisat2" --version | head -1)"
  echo "samtools_version=$("$BIN/samtools" --version | head -1)"
} > "$OUT/run_context.txt"

# --------------------------------------------------- 1. align to target refs
"$BIN/hisat2" -p "$THREADS" -x "$VIRAL_IDX" "${READ_ARGS[@]}" \
  --no-spliced-alignment --no-unal \
  --summary-file "$OUT/viral_alignment.summary.txt" \
  -S "$OUT/viral_hits.sam" 2> "$OUT/hisat2_viral.log"

# --------------------------------------------- 2. collapse unique sequences
"$BIN/python" -I "$SCRIPTS/score_library.py" candidates \
  --viral-sam "$OUT/viral_hits.sam" \
  --out-fasta "$OUT/candidates.fasta" \
  --out-tsv "$OUT/candidates.tsv" 2>> "$OUT/hisat2_viral.log"

# ------------------------------------------------- 3. human homology re-check
if [ -s "$OUT/candidates.fasta" ]; then
  "$BIN/hisat2" -p "$THREADS" -f -x "$HUMAN_IDX" -U "$OUT/candidates.fasta" \
    -k 1 --summary-file "$OUT/human_recheck.summary.txt" \
    -S "$OUT/candidates_vs_human.sam" 2> "$OUT/hisat2_human.log"
else
  : > "$OUT/candidates_vs_human.sam"
fi

# ------------------------------------------------------ 4. apply the rule
"$BIN/python" -I "$SCRIPTS/score_library.py" score \
  --library "$SRR" --cohort "$COHORT" \
  --candidate-table "$OUT/candidates.tsv" \
  --human-sam "$OUT/candidates_vs_human.sam" \
  --library-summary "$OUT/viral_alignment.summary.txt" \
  --accession-table "$ACC_TABLE" \
  --mask-report "$MASK_REPORT" \
  --out-json "$OUT/screen_result.json" \
  --out-tsv "$OUT/screen_result.tsv"

# Keep the target hits as a compressed BAM; drop the large SAM.
"$BIN/samtools" view -b -o "$OUT/viral_hits.bam" "$OUT/viral_hits.sam"
rm -f "$OUT/viral_hits.sam"
gzip -f "$OUT/candidates_vs_human.sam"

echo "finished_utc=$(date -u +%FT%TZ)" >> "$OUT/run_context.txt"
cat "$OUT/screen_result.tsv"
