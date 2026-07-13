#!/bin/bash
# 01_download_gse281367.sbatch  --  Layer A A2 (download arm)
# Download Zhu/Hu 2026 GSE281367 snATAC (PRJNA1183734, 12 donors) from SRA and
# rename SRA reads to 10x cellranger-atac convention.
#
# SRA read structure (verified on SRR31281574, asserted per-run below):
#   _1 = 8 bp  sample index  -> I1
#   _2 = 50 bp genomic read1 -> R1
#   _3 = 16 bp cell barcode  -> R2   (the read ENA's 2-file view drops)
#   _4 = 50 bp genomic read2 -> R3
# cellranger-atac count expects R1 (genomic), R2 (16bp barcode), R3 (genomic).
#
#SBATCH --job-name=fasterqdump
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --array=1-12
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/dl_%A_%a.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/dl_%A_%a.err

set -eo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
module load SRA-Toolkit/3.2.1

META="$PROJ/data/GSE281367/metadata/donor_pairing.csv"
# array task i -> data row (skip header)
ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$META")
DONOR=$(echo "$ROW" | cut -d, -f1)
SRR=$(echo   "$ROW" | cut -d, -f2)
COND=$(echo  "$ROW" | cut -d, -f5)
echo "[$(date)] task=$SLURM_ARRAY_TASK_ID donor=$DONOR srr=$SRR cond=$COND host=$(hostname)"

OUT="$PROJ/data/raw/GSE281367/fastq/$DONOR"
CACHE="$PROJ/data/raw/GSE281367/sra_cache"
TMP="$PROJ/data/raw/GSE281367/tmp_${SRR}"
mkdir -p "$OUT" "$CACHE" "$TMP"

# skip if already finalized
if [ -f "$OUT/${DONOR}_S1_L001_R2_001.fastq.gz" ]; then
  echo "[$(date)] $DONOR already downloaded; skipping."; exit 0
fi

# 1) prefetch the .sra (robust for large files); cache on the data volume, not ~/ncbi
prefetch "$SRR" -O "$CACHE" --max-size u
SRA_PATH="$CACHE/$SRR/$SRR.sra"; [ -f "$SRA_PATH" ] || SRA_PATH="$CACHE/$SRR.sra"

# 2) extract ALL reads (biological + technical barcode/index)
fasterq-dump "$SRA_PATH" --split-files --include-technical -e "$SLURM_CPUS_PER_TASK" \
  -t "$TMP" -O "$OUT"

# 3) map reads to cellranger-atac names BY READ LENGTH (robust to SRA read ordering),
#    asserting the expected {8,16,50,50} structure so any deviant run fails loudly.
declare -A ROLE
GEN=()   # the two ~50bp genomic reads, kept in SRA file order for R1/R3
for f in "$OUT/${SRR}"_*.fastq; do
  L=$(awk 'NR==2{print length($0); exit}' "$f")
  if   [ "$L" -eq 16 ]; then ROLE["$f"]=R2
  elif [ "$L" -le 12 ]; then ROLE["$f"]=I1          # 8bp sample index
  elif [ "$L" -ge 40 ]; then GEN+=("$f")            # genomic
  else echo "ERROR: unexpected read length $L in $f" >&2; exit 3
  fi
done
if [ "${#GEN[@]}" -ne 2 ]; then
  echo "ERROR: expected 2 genomic reads, found ${#GEN[@]} for $SRR" >&2; exit 4
fi
# GEN sorted by SRA filename order -> R1 (first) then R3 (second)
IFS=$'\n' GEN=($(printf '%s\n' "${GEN[@]}" | sort)); unset IFS
ROLE["${GEN[0]}"]=R1; ROLE["${GEN[1]}"]=R3

have_r1=0; have_r2=0; have_r3=0
for f in "${!ROLE[@]}"; do
  R=${ROLE[$f]}
  mv "$f" "$OUT/${DONOR}_S1_L001_${R}_001.fastq"
  [ "$R" = R1 ] && have_r1=1; [ "$R" = R2 ] && have_r2=1; [ "$R" = R3 ] && have_r3=1
done
if [ "$have_r1$have_r2$have_r3" != "111" ]; then
  echo "ERROR: missing one of R1/R2/R3 for $SRR (r1=$have_r1 r2=$have_r2 r3=$have_r3)" >&2; exit 5
fi

# 4) compress (parallel pigz from the rnaseq env; else pigz on PATH; else gzip); clean up .sra + temp
PIGZBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pigz
if [ -x "$PIGZBIN" ]; then GZ="$PIGZBIN -p $SLURM_CPUS_PER_TASK"
elif command -v pigz >/dev/null; then GZ="pigz -p $SLURM_CPUS_PER_TASK"
else GZ="gzip"; fi
echo "[$(date)] compressing with: $GZ"
$GZ "$OUT"/${DONOR}_S1_L001_*.fastq
rm -rf "$TMP" "$CACHE/$SRR" "$SRA_PATH" 2>/dev/null || true

echo "[$(date)] $DONOR done. Files:"
ls -la "$OUT"
