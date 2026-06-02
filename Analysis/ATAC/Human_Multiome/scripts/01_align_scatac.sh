#!/bin/bash
#SBATCH --job-name=scatac_align
#SBATCH --array=1-18
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=36:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/align_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome/logs/align_%A_%a.err

# =============================================================================
# Module 2 Step 1: scATAC-seq Alignment + Fragment Generation
# GSE244832 uses combinatorial indexing (NOT 10x Chromium)
# Barcodes are embedded in FASTQ headers by ATACdemultiplex
# Pipeline: preprocess FASTQs → bowtie2 align → filter → generate fragments
# =============================================================================

set -euo pipefail

# ---- Paths ----
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_DIR="$PROJECT_ROOT/Analysis/ATAC/Human_Multiome"
METADATA="$PROJECT_ROOT/data/GSE244832/metadata/donor_pairing.csv"
FQDIR="$PROJECT_ROOT/data/GSE244832/fastq"
OUTDIR="$PIPELINE_DIR/results"
BT2_INDEX="/gpfs/commons/home/jameslee/reference_genome/bowtie2/GRCh38/GRCh38"
BLACKLIST="/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed"
THREADS=${SLURM_CPUS_PER_TASK:-16}

# ---- Parse sample from array task ID ----
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$METADATA")
DONOR_ID=$(echo "$LINE" | cut -d, -f1)
SRR=$(echo "$LINE" | cut -d, -f3)
CONDITION=$(echo "$LINE" | cut -d, -f8)

# Auto-detect FASTQ naming convention (10x-style vs SRA-style)
if [[ -f "$FQDIR/$SRR/${SRR}_S1_L001_R1_001.fastq.gz" ]]; then
    R1="$FQDIR/$SRR/${SRR}_S1_L001_R1_001.fastq.gz"
    R2="$FQDIR/$SRR/${SRR}_S1_L001_R2_001.fastq.gz"
elif [[ -f "$FQDIR/$SRR/${SRR}_1.fastq.gz" ]]; then
    R1="$FQDIR/$SRR/${SRR}_1.fastq.gz"
    R2="$FQDIR/$SRR/${SRR}_2.fastq.gz"
else
    echo "ERROR: No FASTQ files found for $SRR in $FQDIR/$SRR/"
    ls "$FQDIR/$SRR/" 2>/dev/null
    exit 1
fi

echo "============================================"
echo "scATAC Alignment: $DONOR_ID ($SRR) [$CONDITION]"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "CPUs: $THREADS, Mem: ${SLURM_MEM_PER_NODE:-96G}"
echo "============================================"

# Validate inputs
for f in "$R1" "$R2" "$BLACKLIST"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing file: $f"
        exit 1
    fi
done

# ---- Load modules ----
module purge 2>/dev/null
module load Bowtie2/2.5.4-linux-x86_64
module load SAMtools/1.21
module load BEDTools/2.31.0-GCC-12.3.0

# ---- Create output dirs ----
mkdir -p "$OUTDIR/alignment" "$OUTDIR/fragments" "$OUTDIR/qc"
TMPDIR_SAMPLE="$OUTDIR/alignment/${DONOR_ID}_tmp"
mkdir -p "$TMPDIR_SAMPLE"

# ---- Step 1: Align with bowtie2 ----
# Barcode is embedded in FASTQ header comment: @SRR.N BARCODE:instrument:...
# Use process substitution to prepend barcode to read name on the fly
# Output read name: BARCODE:SRR.N (barcode extracted from before first colon in comment)
echo ""
echo "[$(date +%H:%M:%S)] Step 1: Aligning with bowtie2 (barcode-aware)..."

bowtie2 \
    --very-sensitive \
    -X 2000 \
    --no-mixed \
    --no-discordant \
    -p "$THREADS" \
    -x "$BT2_INDEX" \
    -1 <(zcat "$R1" | awk 'NR%4==1 { split($2,bc,":"); sub(/^@/,"@"bc[1]":") } {print}') \
    -2 <(zcat "$R2" | awk 'NR%4==1 { split($2,bc,":"); sub(/^@/,"@"bc[1]":") } {print}') \
    2> "$OUTDIR/qc/${DONOR_ID}_bowtie2.log" \
| samtools sort -@ 4 -m 4G -T "$TMPDIR_SAMPLE/sort" \
    -o "$OUTDIR/alignment/${DONOR_ID}.sorted.bam"

samtools index -@ 4 "$OUTDIR/alignment/${DONOR_ID}.sorted.bam"

echo "[$(date +%H:%M:%S)] Alignment complete."
grep "overall alignment rate" "$OUTDIR/qc/${DONOR_ID}_bowtie2.log" || true

# ---- Step 2: Filter BAM ----
# Keep: proper pairs (-f 0x2), MAPQ>=30 (-q 30)
# Remove: unmapped, secondary, failed QC, supplementary (-F 0xB0C)
# Remove: chrM reads, ENCODE blacklist regions
echo ""
echo "[$(date +%H:%M:%S)] Step 2: Filtering BAM..."

# Get all chromosomes except chrM and unplaced
CHROMS=$(samtools idxstats "$OUTDIR/alignment/${DONOR_ID}.sorted.bam" \
    | awk '$1 != "chrM" && $1 != "*" && $3 > 0 {printf $1" "}')

samtools view -b -f 0x2 -F 0xB0C -q 30 \
    "$OUTDIR/alignment/${DONOR_ID}.sorted.bam" $CHROMS \
| bedtools intersect -v -abam stdin -b "$BLACKLIST" \
> "$OUTDIR/alignment/${DONOR_ID}.filtered.bam"

samtools index -@ 4 "$OUTDIR/alignment/${DONOR_ID}.filtered.bam"

echo "[$(date +%H:%M:%S)] Filtering complete."

# ---- Step 3: Generate fragment file ----
# Extract fragment coordinates from properly paired reads
# Fragment = (chr, start, end, barcode, count)
# Only process reads where TLEN > 0 (avoids double-counting pairs)
# Collapse identical fragments (same position + barcode) and count
echo ""
echo "[$(date +%H:%M:%S)] Step 3: Generating fragment file..."

samtools view "$OUTDIR/alignment/${DONOR_ID}.filtered.bam" \
| awk '$9 > 0 {
    n = index($1, ":");
    barcode = substr($1, 1, n-1);
    chr = $3;
    start = $4 - 1;
    end = start + $9;
    if (end > start) print chr "\t" start "\t" end "\t" barcode
}' \
| sort -k1,1 -k2,2n -k3,3n -k4,4 -S 16G --parallel=4 -T "$TMPDIR_SAMPLE" \
| awk '{
    key = $1"\t"$2"\t"$3"\t"$4;
    if (key == prev) { count++ }
    else { if (NR > 1) print prev "\t" count; prev = key; count = 1 }
} END { if (NR > 0) print prev "\t" count }' \
| bgzip > "$OUTDIR/fragments/${DONOR_ID}_fragments.tsv.gz"

tabix -p bed "$OUTDIR/fragments/${DONOR_ID}_fragments.tsv.gz"

echo "[$(date +%H:%M:%S)] Fragment file generated."

# ---- Step 4: QC statistics ----
echo ""
echo "[$(date +%H:%M:%S)] Step 4: Computing QC statistics..."

TOTAL_READS=$(samtools view -c "$OUTDIR/alignment/${DONOR_ID}.sorted.bam")
MAPPED_READS=$(samtools view -c -F 4 "$OUTDIR/alignment/${DONOR_ID}.sorted.bam")
PROPER_PAIRS=$(samtools view -c -f 2 "$OUTDIR/alignment/${DONOR_ID}.sorted.bam")
FILTERED_READS=$(samtools view -c "$OUTDIR/alignment/${DONOR_ID}.filtered.bam")
TOTAL_FRAGS=$(zcat "$OUTDIR/fragments/${DONOR_ID}_fragments.tsv.gz" | wc -l)
UNIQUE_FRAGS=$TOTAL_FRAGS  # already collapsed
N_BARCODES=$(zcat "$OUTDIR/fragments/${DONOR_ID}_fragments.tsv.gz" | cut -f4 | sort -u -S 4G | wc -l)
MEDIAN_FRAGS=$(zcat "$OUTDIR/fragments/${DONOR_ID}_fragments.tsv.gz" \
    | cut -f4 | sort | uniq -c | awk '{print $1}' | sort -n \
    | awk '{a[NR]=$1} END {if (NR%2) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2}')

# Write QC CSV
echo "donor_id,srr,condition,total_reads,mapped_reads,proper_pairs,filtered_reads,unique_fragments,n_barcodes,median_frags_per_barcode" \
    > "$OUTDIR/qc/${DONOR_ID}_stats.csv"
echo "${DONOR_ID},${SRR},${CONDITION},${TOTAL_READS},${MAPPED_READS},${PROPER_PAIRS},${FILTERED_READS},${UNIQUE_FRAGS},${N_BARCODES},${MEDIAN_FRAGS}" \
    >> "$OUTDIR/qc/${DONOR_ID}_stats.csv"

echo ""
echo "  Total reads:     $TOTAL_READS"
echo "  Mapped reads:    $MAPPED_READS"
echo "  Proper pairs:    $PROPER_PAIRS"
echo "  After filtering: $FILTERED_READS"
echo "  Unique fragments: $UNIQUE_FRAGS"
echo "  Barcodes (cells): $N_BARCODES"
echo "  Median frags/cell: $MEDIAN_FRAGS"

# ---- Cleanup ----
# Keep filtered BAM and fragments; remove raw sorted BAM to save space
rm -f "$OUTDIR/alignment/${DONOR_ID}.sorted.bam" "$OUTDIR/alignment/${DONOR_ID}.sorted.bam.bai"
rm -rf "$TMPDIR_SAMPLE"

echo ""
echo "============================================"
echo "Completed: $DONOR_ID ($SRR) - $(date)"
echo "Outputs:"
echo "  BAM: $OUTDIR/alignment/${DONOR_ID}.filtered.bam"
echo "  Fragments: $OUTDIR/fragments/${DONOR_ID}_fragments.tsv.gz"
echo "  QC: $OUTDIR/qc/${DONOR_ID}_stats.csv"
echo "============================================"
