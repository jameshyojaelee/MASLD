#!/bin/bash
#SBATCH --job-name=h25_stage1
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# H25 stage 1: per run, download (own directory), bait-select read pairs near the 95 target
# SNVs, align to the whole genome, remap filter, random dedup, fragment allele counts.
# Usage: sbatch --output=<log> run_stage1.sbatch <run_list> <data_dir> <targets_dir>
set -euo pipefail

RUNLIST=$1
DATA=$2
TDIR=$3
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H25_cutrun_allelic_imbalance
PY="/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/python -I"
IDX=/gpfs/commons/home/jameslee/reference_genome/bowtie2/GRCh38/GRCh38
BBDUK=/gpfs/commons/home/jameslee/bbmap/bbduk.sh
BT2OPT="--end-to-end --very-sensitive -X 1000 --no-mixed --no-discordant --no-unal"
T=${SLURM_CPUS_PER_TASK:-16}

module load Bowtie2/2.5.4-linux-x86_64 SAMtools/1.18-GCC-12.3.0 SRA-Toolkit/3.2.1 Java/17

echo "host $(hostname) start $(date -u)"
for SRR in $(cat "$RUNLIST"); do
  OUT=$DATA/stage1/$SRR
  if [[ -f $OUT/DONE ]]; then echo "$SRR done, skip"; continue; fi
  mkdir -p "$OUT" "$DATA/sra/$SRR" "$DATA/tmp/$SRR"
  SRA=$DATA/sra/$SRR/$SRR/$SRR.sra
  if [[ ! -s $SRA ]]; then
    for i in 1 2 3 4 5; do
      prefetch "$SRR" -O "$DATA/sra/$SRR" --max-size 20G && break
      echo "prefetch retry $i"; sleep 60
    done
  fi
  vdb-validate "$SRA" > "$OUT/vdb_validate.txt" 2>&1
  fasterq-dump "$SRA" --split-files -e "$T" -t "$DATA/tmp/$SRR" -O "$DATA/tmp/$SRR"
  head -8 "$DATA/tmp/$SRR/${SRR}_1.fastq" > "$OUT/fastq_head.txt"
  "$BBDUK" -Xmx8g in1="$DATA/tmp/$SRR/${SRR}_1.fastq" in2="$DATA/tmp/$SRR/${SRR}_2.fastq" \
    ref="$TDIR/bait_windows.fa" k=25 hdist=1 threads="$T" \
    outm1="$OUT/bait_1.fq.gz" outm2="$OUT/bait_2.fq.gz" stats="$OUT/bbduk_stats.txt" 2> "$OUT/bbduk.log"
  rm -f "$DATA/tmp/$SRR/${SRR}_1.fastq" "$DATA/tmp/$SRR/${SRR}_2.fastq"

  bowtie2 $BT2OPT -p "$T" -x "$IDX" -1 "$OUT/bait_1.fq.gz" -2 "$OUT/bait_2.fq.gz" 2> "$OUT/bowtie2.log" \
    | samtools view -b -f 2 -F 0x904 -q 30 - \
    | samtools sort -n -@ 4 -m 1G -o "$OUT/filt.nsort.bam" -
  $PY "$S/wasp_remap.py" make --bam "$OUT/filt.nsort.bam" --sites "$TDIR/targets.tsv" --out-prefix "$OUT/wasp"
  bowtie2 $BT2OPT -p "$T" -x "$IDX" -1 "$OUT/wasp.remap_1.fq" -2 "$OUT/wasp.remap_2.fq" 2> "$OUT/bowtie2_remap.log" \
    | samtools view -b -o "$OUT/wasp.remapped.bam" -
  $PY "$S/wasp_remap.py" check --bam "$OUT/filt.nsort.bam" --out-prefix "$OUT/wasp" \
    --remapped "$OUT/wasp.remapped.bam" --out-bam "$OUT/kept.nsort.bam"
  samtools sort -@ 4 -m 1G -o "$OUT/kept.bam" "$OUT/kept.nsort.bam"
  samtools index "$OUT/kept.bam"
  $PY "$S/count_sites.py" --bam "$OUT/kept.bam" --sites "$TDIR/targets.tsv" --seed 20260923 \
    --sample "$SRR" --out "$OUT/counts.tsv"
  rm -f "$OUT/wasp.remap_1.fq" "$OUT/wasp.remap_2.fq" "$OUT/kept.nsort.bam"
  touch "$OUT/DONE"
  echo "$SRR finished $(date -u)"
done
echo "end $(date -u)"
