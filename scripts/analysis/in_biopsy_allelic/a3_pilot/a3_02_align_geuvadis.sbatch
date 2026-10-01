#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --array=0-29
# A3 pilot: align one Geuvadis RNA-seq library (30 unrelated EUR founders) to the
# primary-chromosome STAR index, two-pass, with a read group. Tasks 0-4 also run
# STAR WASP mode with the individual's 1kGP DNA het SNVs, to check that
# --waspOutputMode works together with --twopassMode Basic (PRESPEC_A gate G1).
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out a3_02_align_geuvadis.sbatch
set -euo pipefail
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external
idx=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/star_2.7.11b
env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1

line=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR == i + 2' "$ext/geuvadis/selected_30_eur.tsv")
run=$(printf '%s\n' "$line" | cut -f1); ind=$(printf '%s\n' "$line" | cut -f2)
fq=$ext/geuvadis/fastq; out=$ext/geuvadis/bam/$run
mkdir -p "$out"
test -s "$fq/${run}_1.fastq.gz" && test -s "$fq/${run}_2.fastq.gz"

wasp=()
if [ "$SLURM_ARRAY_TASK_ID" -lt 5 ]; then
  panel=$ext/allelic_refs/kgp_highcov
  for c in $(seq 1 22); do
    bcftools view -s "$ind" -v snps -m2 -M2 -g het -Oz -o "$out/tmp.chr$c.vcf.gz" \
      "$panel/1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  done
  # shellcheck disable=SC2046
  bcftools concat -Ov -o "$out/$ind.dna_het_snvs.vcf" $(printf "$out/tmp.chr%s.vcf.gz " $(seq 1 22))
  rm "$out"/tmp.chr*.vcf.gz
  wasp=(--waspOutputMode SAMtag --varVCFfile "$out/$ind.dna_het_snvs.vcf")
fi

"$env/STAR" --runThreadN 12 --genomeDir "$idx" --twopassMode Basic \
  --readFilesIn "$fq/${run}_1.fastq.gz" "$fq/${run}_2.fastq.gz" --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
  --outSAMattributes NH HI AS nM MD --outSAMattrRGline "ID:$run" "SM:$ind" \
  --outFileNamePrefix "$out/$run." "${wasp[@]}"
"$env/samtools" index "$out/$run.Aligned.sortedByCoord.out.bam"
if [ "${#wasp[@]}" -gt 0 ]; then
  "$env/samtools" view "$out/$run.Aligned.sortedByCoord.out.bam" | awk '
    {for (i = 12; i <= NF; i++) if ($i ~ /^vW:i:/) {t[$i]++; n++}} END {for (k in t) print k, t[k]; print "reads_with_vW", n + 0}' \
    > "$out/$run.wasp_tag_counts.txt"
fi
echo "$run $ind done"
