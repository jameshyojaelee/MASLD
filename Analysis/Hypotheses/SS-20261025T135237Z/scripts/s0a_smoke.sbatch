#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=interactive
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS job S0a: smoke test of the SteatoSITE alignment on 2 samples (S0005: 2 runs, 176M pairs; S1169: 4 runs).
# STAR flags are the release run_star rule (GSE135251 Snakefile:94-112) except: --runThreadN 8, lanes comma-joined
# with one @RG per run. No trimming (release had none; STAR soft-clips the Pico v2 read-2 adapter bases).
# Measures strandedness (ReadsPerGene cols 3 vs 4), time, peak memory, BAM size. Then, only if read 1 is antisense
# (featureCounts -s 2, the strand the r8 counter enforces), runs the unchanged r8 verified counter.
# Expression values are not inspected here; only strand and technical totals.
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
FQ=$ROOT/data/external/SteatoSite/fastq
META=$FQ/PRJEB58625_run_metadata.tsv
BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
IDX=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/star_index
OUT=$SS/out_s0a_smoke
[ -e "$OUT" ] && { echo "refusing to overwrite $OUT"; exit 1; }
mkdir -p "$OUT"
sha256sum "$0" > "$OUT/script_sha256.txt" 2>/dev/null || true
for S in S0005 S1169; do
  d=$OUT/$S; mkdir -p "$d"
  runs=$(awk -F'\t' -v s=$S 'NR>1 && $5==s {print $1}' "$META" | sort)
  r1=$(for r in $runs; do printf "%s," "$FQ/${r}_1.fastq.gz"; done | sed 's/,$//')
  r2=$(for r in $runs; do printf "%s," "$FQ/${r}_2.fastq.gz"; done | sed 's/,$//')
  rg=$(for r in $runs; do printf "ID:%s SM:%s , " "$r" "$S"; done | sed 's/ , $//')
  echo "$S runs: $runs"
  /usr/bin/time -v -o "$d/time_star.txt" $BIN/STAR --runThreadN 8 --genomeDir $IDX \
    --readFilesIn "$r1" "$r2" --readFilesCommand zcat --outSAMattrRGline $rg \
    --outFileNamePrefix "$d/$S." --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic --quantMode GeneCounts --limitBAMsortRAM 30000000000 > "$d/star.log" 2>&1
  $BIN/samtools index -@ 8 "$d/$S.Aligned.sortedByCoord.out.bam"
  ls -l "$d/$S.Aligned.sortedByCoord.out.bam" >> "$OUT/bam_sizes.txt"
  # strand: rows 1-4 of ReadsPerGene are N_* summaries; col2 unstranded, col3 read1 sense, col4 read1 antisense
  awk -v s=$S 'NR>4{a+=$3; b+=$4} END{printf "%s\tcol3_read1_sense=%d\tcol4_read1_antisense=%d\tfrac_col4=%.4f\n", s, a, b, b/(a+b)}' \
    "$d/$S.ReadsPerGene.out.tab" >> "$OUT/strandedness.tsv"
  grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads unmapped: too short|Mismatch rate" "$d/$S.Log.final.out" | sed "s/^/$S /" >> "$OUT/star_summary.txt"
done
cat "$OUT/strandedness.tsv"
ok=$(awk -F'\t' '{split($4,a,"="); if (a[2] > 0.8) n++} END{print (n==2) ? "yes" : "no"}' "$OUT/strandedness.tsv")
echo "read1 antisense in both (>0.8): $ok" | tee "$OUT/strand_gate.txt"
if [ "$ok" == yes ]; then
  for S in S0005 S1169; do
    export SOURCE_BAM="$OUT/$S/$S.Aligned.sortedByCoord.out.bam" SAMPLE_ID=$S LIBRARY_LAYOUT=PAIRED \
      ANNOTATION_ROOT=$ROOT/figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r8 \
      COUNT_OUTPUT_ROOT="$OUT/$S/verified_count"
    unset SLURM_TMPDIR
    export COUNT_SCRATCH_ROOT=$(mktemp -d "$OUT/$S/count_scratch-XXXXXX")
    /usr/bin/time -v -o "$OUT/$S/time_count.txt" bash $ROOT/scripts/analysis/alt_contig_recount/run_verified_count.sbatch
    ls "$COUNT_OUTPUT_ROOT"
  done
fi
echo "S0a done"
