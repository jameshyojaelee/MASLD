#!/bin/bash
#SBATCH --job-name=recount_vm38
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/recount_vm38_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/recount_vm38_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00

# Phase A: Re-quantify InHouse_MCD + Public_MCD BAMs with GENCODE vM38 GTF
#
# Why this works: All samples were aligned to GRCm39 genome assembly.
# Only the GTF annotation differs (vM33 vs vM38). featureCounts counts
# against the GTF, not the genome, so re-running featureCounts on existing
# BAMs with the vM38 GTF gives correct vM38 counts without re-alignment.
#
# Strandedness:
#   InHouse_MCD:  reverse-stranded (-s 2)  [TruSeq Stranded]
#   GSE156918:    unstranded (-s 0)
#   GSE205974:    unstranded (-s 0)

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

MOUSE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse"
OUTDIR="$MOUSE/Unified_Integration"
GTF="$MOUSE/Public_Diet_Models/reference/raw/gencode.vM38.annotation.gtf"
THREADS=16

mkdir -p "$OUTDIR/counts/featurecounts" "$OUTDIR/logs"

echo "=== Phase A: featureCounts Re-Quantification with vM38 GTF ==="
echo "GTF: $GTF"
echo "Start: $(date)"
echo ""

# --- 1) InHouse MCD: paired-end, reverse-stranded (-s 2) ---
echo "=== 1/3: InHouse MCD (PE, -s 2) ==="
BAMS_INHOUSE=()
for sample in 1CF 1CM 1MF 1MM 2CF 2CM 2MF 2MM 3CF 3CM 3MF 3MM; do
  bam="$MOUSE/InHouse_MCD/alignments/star/${sample}/${sample}.Aligned.sortedByCoord.out.bam"
  if [ -f "$bam" ]; then
    BAMS_INHOUSE+=("$bam")
  else
    echo "  WARNING: Missing BAM: $bam"
  fi
done
echo "  BAMs found: ${#BAMS_INHOUSE[@]}/12"

featureCounts \
  -T "$THREADS" \
  -p --countReadPairs \
  -s 2 \
  -a "$GTF" \
  -o "$OUTDIR/counts/featurecounts/gene_counts_inhouse.txt" \
  "${BAMS_INHOUSE[@]}" 2>&1

echo "  Done"
echo ""

# --- 2) GSE156918: single-end, unstranded (-s 0) ---
echo "=== 2/3: GSE156918 (SE, -s 0) ==="
BAMS_GSE156918=()
for bam in "$MOUSE"/Public_MCD/GSE156918/alignments/star/GSM*/GSM*.Aligned.sortedByCoord.out.bam; do
  if [ -f "$bam" ]; then
    BAMS_GSE156918+=("$bam")
  fi
done
echo "  BAMs found: ${#BAMS_GSE156918[@]}"

featureCounts \
  -T "$THREADS" \
  -s 0 \
  -a "$GTF" \
  -o "$OUTDIR/counts/featurecounts/gene_counts_gse156918.txt" \
  "${BAMS_GSE156918[@]}" 2>&1

echo "  Done"
echo ""

# --- 3) GSE205974: paired-end, unstranded (-s 0) ---
echo "=== 3/3: GSE205974 (PE, -s 0) ==="
BAMS_GSE205974=()
for bam in "$MOUSE"/Public_MCD/GSE205974/alignments/star/GSM*/GSM*.Aligned.sortedByCoord.out.bam; do
  if [ -f "$bam" ]; then
    BAMS_GSE205974+=("$bam")
  fi
done
echo "  BAMs found: ${#BAMS_GSE205974[@]}"

featureCounts \
  -T "$THREADS" \
  -p --countReadPairs \
  -s 0 \
  -a "$GTF" \
  -o "$OUTDIR/counts/featurecounts/gene_counts_gse205974.txt" \
  "${BAMS_GSE205974[@]}" 2>&1

echo "  Done"
echo ""

TOTAL=$(( ${#BAMS_INHOUSE[@]} + ${#BAMS_GSE156918[@]} + ${#BAMS_GSE205974[@]} ))
echo "=== Phase A complete: $(date) ==="
echo "Total samples re-quantified: $TOTAL"
echo "Outputs:"
echo "  InHouse:   $OUTDIR/counts/featurecounts/gene_counts_inhouse.txt"
echo "  GSE156918: $OUTDIR/counts/featurecounts/gene_counts_gse156918.txt"
echo "  GSE205974: $OUTDIR/counts/featurecounts/gene_counts_gse205974.txt"
