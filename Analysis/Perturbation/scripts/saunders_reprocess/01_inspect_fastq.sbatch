#!/bin/bash
#SBATCH --job-name=prefetch
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Perturbation/scripts/saunders_reprocess/logs/01_inspect_%j.log

# Step 1: Download smallest SRR from PRJNA1151062 (Saunders 2025 Perturb-Multimodal)
# and characterize FASTQ structure to determine chemistry layout.
#
# Target SRR30344305 (676 MB lite, GSM8478327, ~5M reads of 202 avg-length = paired ~101bp)
# Goal: decide between cellranger multi (Flex Feature Barcode Kit) vs custom demux.

set -euo pipefail

module load SRA-Toolkit/3.2.1

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRA_DIR=$PROJ/data/perturbation/datasets/saunders2025/raw/sra
FQ_DIR=$PROJ/data/perturbation/datasets/saunders2025/raw/fastq
SCRATCH=/scratch/jameslee/saunders_sra
mkdir -p "$SRA_DIR" "$FQ_DIR" "$SCRATCH"

ACCESSIONS=(SRR30344305 SRR30344303)   # 2 smallest; one is enough for inspection

cd "$SCRATCH"

for ACC in "${ACCESSIONS[@]}"; do
    echo "=== $ACC ==="
    echo "[$(date +%T)] prefetch $ACC"
    prefetch --output-directory "$SRA_DIR" --max-size 50g "$ACC"

    SRA_FILE=$(ls "$SRA_DIR/$ACC"/*.sra 2>/dev/null | head -1 || true)
    if [ -z "$SRA_FILE" ]; then
        echo "ERROR: no .sra file for $ACC"; exit 1
    fi
    echo "[$(date +%T)] fasterq-dump $ACC ($(du -sh "$SRA_FILE"))"

    fasterq-dump \
        --split-files \
        --include-technical \
        --threads "${SLURM_CPUS_PER_TASK:-8}" \
        --temp "$SCRATCH/tmp_$ACC" \
        --outdir "$FQ_DIR" \
        "$SRA_FILE"

    echo "[$(date +%T)] FASTQ files for $ACC:"
    ls -la "$FQ_DIR"/${ACC}_* 2>/dev/null || true
done

# Inspection summary per SRR
for ACC in "${ACCESSIONS[@]}"; do
    echo ""
    echo "########## $ACC inspection ##########"
    for FQ in "$FQ_DIR"/${ACC}_*.fastq; do
        [ -e "$FQ" ] || continue
        BN=$(basename "$FQ")
        N_LINES=$(wc -l < "$FQ")
        N_READS=$((N_LINES / 4))
        # Read length distribution (first 10k reads to keep it fast)
        LEN_DIST=$(awk 'NR%4==2{print length($0)}' "$FQ" | head -10000 | sort -n | uniq -c | awk '{printf "%s=%s  ", $2, $1}')
        echo "--- $BN ($N_READS reads) ---"
        echo "Read-length distribution (first 10k): $LEN_DIST"
        echo "First 4 reads:"
        head -16 "$FQ"
        # Look for known sgRNA constant-region motif (mmc6 probe scaffold)
        SEEN_CONST=$(head -2000000 "$FQ" | awk 'NR%4==2' | grep -c "GCTATGCTGTTTCCAGCTTAGCTCT" || true)
        echo "Reads (of first 500k) containing sgRNA-constant-region GCTATGCTGTTTCCAGCTTAGCTCT: $SEEN_CONST"
        echo ""
    done

    # gzip the FASTQs for downstream use
    echo "[$(date +%T)] gzip ${ACC} fastq files"
    gzip -f "$FQ_DIR"/${ACC}_*.fastq || true
done

echo "=== DONE at $(date) ==="
