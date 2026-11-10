#!/bin/bash
# SS S1 per-library unit: STAR (release run_star flags + --runThreadN 8, one @RG per lane run), samtools index, then the
# unchanged r8 verified fragment counter (-s 2). Mates are put in kit order per run first (PRESPEC_SS amendment A1:
# ENA read labels are swapped in some runs; per-run STAR probe, s1_probe.sh), and counting runs only if the
# aligned library has read 1 antisense in > 0.80 of stranded gene reads (STAR ReadsPerGene). Called by the nslab array and by
# the interactive workers on the same sample list; a per-sample lock (mkdir) prevents two jobs working on one library,
# and a lock whose owning job is no longer in the queue is treated as stale. Usage: s1_one.sh <sample_alias>
set -euo pipefail
S=$1
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
FQ=$ROOT/data/external/SteatoSite/fastq
META=$FQ/PRJEB58625_run_metadata.tsv
BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
IDX=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/star_index
OUTROOT=$ROOT/RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star
LOCKS=$SS/locks; DONE=$SS/done; mkdir -p "$LOCKS" "$DONE" "$OUTROOT" "$SS/out_s1b_counts/samples"
ME="${SLURM_JOB_ID:-nojob}"
[ -e "$DONE/$S" ] && [ -e "$SS/geno_done/$S" ] && [ -e "$SS/out_s1d_ig/$S/ig_fragments.tsv" ] && { echo "done already: $S"; exit 0; }
if ! mkdir "$LOCKS/$S" 2>/dev/null; then
  owner=$(cat "$LOCKS/$S/owner" 2>/dev/null || echo "")
  if [ -n "$owner" ] && squeue -h -j "$owner" 2>/dev/null | grep -q .; then echo "locked by $owner: $S"; exit 0; fi
  rm -rf "$LOCKS/$S"; mkdir "$LOCKS/$S" || { echo "lock race: $S"; exit 0; }
fi
echo "$ME" > "$LOCKS/$S/owner"
trap 'rm -rf "$LOCKS/$S"' EXIT
d=$OUTROOT/$S
# Library = 11-character code in the submitted FASTQ name (batch code + library number). A sample with two libraries
# (8 samples, two batches each) keeps the library with more read pairs (PRESPEC_SS amendment A2); the other is not aligned.
runs=$(awk -F'\t' -v s=$S 'NR>1 && $5==s {split($22,f,";"); n=split(f[1],p,"/"); match(p[n], /_[0-9][0-9][0-9][0-9][0-9][A-Z][A-Z][0-9][0-9][0-9][0-9]L0/);
  lib=substr(p[n],RSTART+1,11); rc[lib]+=$17; rl[lib]=rl[lib]" "$1}
  END{best=""; for (l in rc) if (best=="" || rc[l]>rc[best] || (rc[l]==rc[best] && l<best)) best=l; print rl[best]}' "$META" | tr ' ' '\n' | grep . | sort)
[ -n "$runs" ] || { echo "no runs for $S"; exit 1; }
# An alignment counts as finished only with mate_order.tsv beside it (written last); anything else is set aside.
if [ -d "$d" ] && ! { [ -s "$d/mate_order.tsv" ] && [ -s "$d/$S.Aligned.sortedByCoord.out.bam.bai" ] && [ -s "$d/$S.Log.final.out" ]; }; then
  # aligned before the mate order was fixed, or incomplete: set aside, redo genotypes and Ig too
  mv "$d" "${d}_mateorder_superseded_$(date +%s)"; rm -f "$SS/geno_done/$S"
  [ -e "$SS/out_s1c_genotypes/vcf/$S.vcf.gz" ] && mv "$SS/out_s1c_genotypes/vcf/$S.vcf.gz" "$SS/out_s1c_genotypes/vcf/$S.vcf.gz.superseded_$(date +%s)"
  [ -e "$SS/out_s1d_ig/$S" ] && mv "$SS/out_s1d_ig/$S" "$SS/out_s1d_ig/${S}_superseded_$(date +%s)"
fi
if [ ! -d "$d" ]; then
  mkdir -p "$d"
  pt=$(mktemp -d "/scratch/ss-probe-${ME}-${S}-XXXXXX")
  bash $SS/scripts/s1_probe.sh "$S" "$pt" $runs > "$d/mate_probe.tsv" || { rm -rf -- "$pt"; mv "$d" "${d}_failed_${ME}_$(date +%s)"; exit 1; }
  cp "$pt/probe.Log.final.out" "$d/mate_probe.Log.final.out"; rm -rf -- "$pt"
  order=$(cut -f1,5 "$d/mate_probe.tsv")
  if [ $(echo "$order" | grep -c .) -ne $(echo "$runs" | wc -w) ] || echo "$order" | grep -q ambiguous; then
    echo "ambiguous mate order, library stopped: $S"; cat "$d/mate_probe.tsv"; mv "$d" "${d}_ambiguous_${ME}_$(date +%s)"; exit 1
  fi
  mate() { awk -F'\t' -v r=$1 '$1==r{print $2}' <(echo "$order"); }
  r1=$(for r in $runs; do [ "$(mate $r)" == swap ] && m=2 || m=1; printf "%s," "$FQ/${r}_${m}.fastq.gz"; done | sed 's/,$//')
  r2=$(for r in $runs; do [ "$(mate $r)" == swap ] && m=1 || m=2; printf "%s," "$FQ/${r}_${m}.fastq.gz"; done | sed 's/,$//')
  rg=$(for r in $runs; do printf "ID:%s SM:%s , " "$r" "$S"; done | sed 's/ , $//')
  echo "$S runs: $runs host: $(hostname) job: $ME"
  T=8; [ "${SLURM_CPUS_PER_TASK:-8}" -lt 8 ] && T=$SLURM_CPUS_PER_TASK  # 4-CPU array tasks; alignments do not depend on thread count
  /usr/bin/time -v -o "$d/time_star.txt" $BIN/STAR --runThreadN $T --genomeDir $IDX \
    --readFilesIn "$r1" "$r2" --readFilesCommand zcat --outSAMattrRGline $rg \
    --outFileNamePrefix "$d/$S." --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic --quantMode GeneCounts --limitBAMsortRAM 30000000000 > "$d/star.log" 2>&1
  rm -rf "$d/$S._STARgenome" "$d/$S._STARpass1" "$d/$S._STARtmp"
  $BIN/samtools index -@ $T "$d/$S.Aligned.sortedByCoord.out.bam"
  echo "$order" > "$d/mate_order.tsv"
fi
asf=$(awk 'NR>4{a+=$3; b+=$4} END{printf "%.4f", b/(a+b)}' "$d/$S.ReadsPerGene.out.tab")
printf '%s\t%s\n' "$S" "$asf" > "$d/read1_antisense_fraction.txt"
if awk -v f=$asf 'BEGIN{exit !(f > 0.8)}'; then
  C=$SS/out_s1b_counts/samples/$S
  if [ ! -s "$C/count_completion.json" ]; then
    [ -e "$C" ] && mv "$C" "${C}_failed_${ME}_$(date +%s)"
    scr=$(mktemp -d "/scratch/ss-count-${ME}-${S}-XXXXXX")
    ( export SAMPLE_ID=$S LIBRARY_LAYOUT=PAIRED SOURCE_BAM=$d/$S.Aligned.sortedByCoord.out.bam \
        ANNOTATION_ROOT=$ROOT/figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r8 \
        COUNT_OUTPUT_ROOT=$C COUNT_SCRATCH_ROOT=$scr SLURM_CPUS_PER_TASK=4; unset SLURM_TMPDIR
      bash $ROOT/scripts/analysis/alt_contig_recount/run_verified_count.sbatch ) || { rm -rf -- "$scr"; exit 1; }
    rm -rf -- "$scr"
  fi
  [ -e "$DONE/$S" ] || { touch "$DONE/$S"; echo "done $S"; }
else
  echo "STRAND FAIL $S: read 1 antisense fraction $asf; not counted (genotypes still written)"; touch "$SS/strand_fail_$S"
fi
# Genotypes at the S0d site list, DIG1 job 62 mpileup/call flags unchanged. Libraries aligned before the list existed are
# filled in by the same call later (s1c pass). Marker: geno_done/<S>.
T=$SS/scripts/s1b_targets.tsv; G=$SS/out_s1c_genotypes
if [ -s "$T" ] && [ ! -e "$SS/geno_done/$S" ]; then
  mkdir -p "$G/vcf" "$G/logs" "$SS/geno_done"
  BCF=/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools
  FA=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa
  $BCF mpileup -f "$FA" -R "$T" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -B -d 100000 "$d/$S.Aligned.sortedByCoord.out.bam" 2> "$G/logs/$S.mpileup.err" \
    | $BCF call -m -Oz -o "$G/vcf/$S.vcf.gz.part" 2> "$G/logs/$S.call.err"
  mv "$G/vcf/$S.vcf.gz.part" "$G/vcf/$S.vcf.gz"
  touch "$SS/geno_done/$S"; echo "genotyped $S"
fi
# IGH constant-region fragments for H38 (scripts/s1d_ig_one.sh; region read through the index, about 1 s per BAM).
[ -e "$SS/out_s1d_ig/$S/ig_fragments.tsv" ] || bash $SS/scripts/s1d_ig_one.sh "$S" "$d/$S.Aligned.sortedByCoord.out.bam" "$SS/out_s1d_ig"
