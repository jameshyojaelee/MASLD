#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=24:00:00
# A5, one library per array task (PRESPEC_A section 6). Realign from FASTQ to the
# primary-chromosome index with STAR WASP mode (this individual's het SNVs from the
# standard imputation, P(het) >= 0.9), then keep only small tables:
#   gene counts (featureCounts, the cohort's original flags, primary-chr GTF);
#   allelic counts (GATK ASEReadCounter, MAPQ 255, WASP-passing vW:i:1 reads) at
#     every site that is het in any leave-gene-out batch (downstream, a gene uses
#     only sites het in its own batch);
#   the 3'/5' coverage ratio (RNA integrity proxy).
# The BAM is deleted at the end (1,466 BAMs would need ~6 TB). Files are removed
# by name; nothing is deleted recursively.
# Usage: sbatch --job-name=<n> --array=1-<N> --output=<exec>/logs/%A_%a.out a5_realign_count.sbatch <libraries.tsv> <a4_restricted_dir> <a5_restricted_dir>
#   libraries.tsv (header): run, individual_id, cohort, fastq1, fastq2 ("-" if single-end), fc_flags
set -euo pipefail
libs=$1; a4=$2; out=$3
IFS=$'\t' read -r run ind cohort fq1 fq2 fcflags < <(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR == i + 1' "$libs")
idx=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49
env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
src=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/a5_allelic
module load bcftools/1.21 htslib/1.21 gatk/4.6.1.0-GCCcore-13.3.0-Java-17 >/dev/null 2>&1
mkdir -p "$out/counts" "$out/ase" "$out/qc" "$out/done"
[ -e "$out/done/$run.done" ] && { echo "$run already done"; exit 0; }
w=$out/work/$run; mkdir -p "$w"

# GATK needs a sequence dictionary next to the FASTA (created once)
if [ ! -s "$idx/genome.primary_chr.dict" ]; then
  "$env/samtools" dict -o "$idx/genome.primary_chr.dict.tmp.$SLURM_ARRAY_TASK_ID" "$idx/genome.primary_chr.fa"
  mv -n "$idx/genome.primary_chr.dict.tmp.$SLURM_ARRAY_TASK_ID" "$idx/genome.primary_chr.dict" 2>/dev/null || true
  rm -f "$idx/genome.primary_chr.dict.tmp.$SLURM_ARRAY_TASK_ID"
fi

# het SNVs: WASP set (standard imputation) and counting set (union over masked batches)
for c in $(seq 1 22); do
  bcftools view -s "$ind" -v snps -i 'GT="het" && FORMAT/GP[0:1] >= 0.9' -Oz -o "$w/wasp.chr$c.vcf.gz" "$a4/imputed/chr$c.standard.vcf.gz"
done
# shellcheck disable=SC2046
bcftools concat -Ov -o "$w/wasp_hets.vcf" $(printf "$w/wasp.chr%s.vcf.gz " $(seq 1 22))
n=0
for f in "$a4"/imputed/chr*.batch*.cis.vcf.gz; do
  n=$((n + 1))
  bcftools view -s "$ind" -v snps -i 'GT="het" && FORMAT/GP[0:1] >= 0.9' -Oz -o "$w/cis.$n.vcf.gz" "$f"
done
# shellcheck disable=SC2046
bcftools concat -Ou $(for k in $(seq 1 "$n"); do printf "$w/cis.%s.vcf.gz " "$k"; done) \
  | bcftools sort -T "$w" -Ou | bcftools norm -d all -Oz -o "$w/count_hets.vcf.gz"
bcftools index -t "$w/count_hets.vcf.gz"

reads=("$fq1"); [ "$fq2" != "-" ] && reads+=("$fq2")
"$env/STAR" --runThreadN 12 --genomeDir "$idx/star_2.7.11b" --twopassMode Basic \
  --readFilesIn "${reads[@]}" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 30000000000 --outSAMattributes NH HI AS nM MD \
  --outSAMattrRGline "ID:$run" "SM:$ind" --waspOutputMode SAMtag --varVCFfile "$w/wasp_hets.vcf" \
  --outFileNamePrefix "$w/$run." --outTmpDir "$w/star_tmp"
bam=$w/$run.Aligned.sortedByCoord.out.bam
"$env/samtools" index "$bam"

# shellcheck disable=SC2086
"$env/featureCounts" -T 12 $fcflags -a "$idx/gencode.v49.primary_chr.gtf" -o "$out/counts/$run.featurecounts.txt" "$bam"
"$env/samtools" view -b -q 255 -d vW:1 -o "$w/wasp_pass.bam" "$bam"
"$env/samtools" index "$w/wasp_pass.bam"
gatk ASEReadCounter -R "$idx/genome.primary_chr.fa" -I "$w/wasp_pass.bam" -V "$w/count_hets.vcf.gz" \
  -mmq 255 -mbq 20 -O "$out/ase/$run.ase.tsv"
"$env/python" "$src/coverage_ratio.py" --bam "$bam" --gtf "$idx/gencode.v49.primary_chr.gtf" \
  --counts "$out/counts/$run.featurecounts.txt" --out "$out/qc/$run.coverage_ratio.json"
cp "$w/$run.Log.final.out" "$out/qc/$run.star_log_final.txt"

# clean up this task's files by name (no recursive deletion)
rm -f "$w"/*.vcf "$w"/*.vcf.gz "$w"/*.tbi "$w"/*.bam "$w"/*.bai "$w/$run".Log.* "$w/$run".SJ.out.tab
for d in "$w/$run._STARgenome" "$w/$run._STARpass1" "$w/star_tmp"; do
  [ -d "$d" ] && { rm -f "$d"/* 2>/dev/null || true; rmdir "$d" 2>/dev/null || true; }
done
rmdir "$w" 2>/dev/null || echo "note: $w not empty after cleanup"
touch "$out/done/$run.done"
echo "$run done"
