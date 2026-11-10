#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-657%40
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/star/%a_%A.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/star/%a_%A.err
# SS job S1: STAR per SteatoSITE library (one task per sample alias in scripts/s1_samples.txt). Flags as S0a (release
# run_star rule plus --runThreadN 8 and one @RG per lane run). A task refuses an existing complete output; an
# incomplete output directory is moved aside as <dir>_failed_<jobid>.
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
FQ=$ROOT/data/external/SteatoSite/fastq
META=$FQ/PRJEB58625_run_metadata.tsv
BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
IDX=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/star_index
OUTROOT=$ROOT/RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SS/scripts/s1_samples.txt)
d=$OUTROOT/$S
if [ -s "$d/$S.Log.final.out" ] && [ -s "$d/$S.Aligned.sortedByCoord.out.bam.bai" ]; then echo "complete: $S"; exit 0; fi
[ -e "$d" ] && mv "$d" "${d}_failed_${SLURM_JOB_ID}"
mkdir -p "$d"
runs=$(awk -F'\t' -v s=$S 'NR>1 && $5==s {print $1}' "$META" | sort)
r1=$(for r in $runs; do printf "%s," "$FQ/${r}_1.fastq.gz"; done | sed 's/,$//')
r2=$(for r in $runs; do printf "%s," "$FQ/${r}_2.fastq.gz"; done | sed 's/,$//')
rg=$(for r in $runs; do printf "ID:%s SM:%s , " "$r" "$S"; done | sed 's/ , $//')
echo "$S runs: $runs host: $(hostname)"
/usr/bin/time -v -o "$d/time_star.txt" $BIN/STAR --runThreadN 8 --genomeDir $IDX \
  --readFilesIn "$r1" "$r2" --readFilesCommand zcat --outSAMattrRGline $rg \
  --outFileNamePrefix "$d/$S." --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic --quantMode GeneCounts --limitBAMsortRAM 30000000000 > "$d/star.log" 2>&1
rm -rf "$d/$S._STARgenome" "$d/$S._STARpass1" "$d/$S._STARtmp"
$BIN/samtools index -@ 8 "$d/$S.Aligned.sortedByCoord.out.bam"
echo "done $S"
