#!/bin/bash
#SBATCH --job-name=V2ctl
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# V2 target screen: build the seeded synthetic positive control and screen it
# with exactly the same script and thresholds as a real library.
#
# Usage: sbatch run_positive_control.sbatch <RESULT_DIR>

set -euo pipefail

RESULT_DIR="${1:?result dir required}"
BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V2_viral_screen
RAW=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw
export PATH="$BIN:$PATH"

HOST_FASTQ="$RAW/GSE162694/fastq/SRR13199872_1.fastq.gz"
CONTIG=LC519824.1          # class-B target reference
CTL_DIR="$RESULT_DIR/positive_control"
mkdir -p "$CTL_DIR"

echo "[$(date -u +%FT%TZ)] building control"
"$BIN/python" -I "$SCRIPTS/make_positive_control.py" \
  --host-fastq "$HOST_FASTQ" \
  --viral-fasta "$RESULT_DIR/downloads_untrusted/viral_refs/viral_refs.fasta" \
  --contig "$CONTIG" \
  --out-fastq "$CTL_DIR/synthetic_positive_control.fastq.gz" \
  --manifest "$CTL_DIR/synthetic_positive_control.manifest.json"

echo "[$(date -u +%FT%TZ)] screening control"
bash "$SCRIPTS/screen_library.sh" \
  "SPIKE_HBV_${CONTIG}" SYNTHETIC_CONTROL "$RESULT_DIR" \
  "$CTL_DIR/synthetic_positive_control.fastq.gz"

echo "[$(date -u +%FT%TZ)] done"
