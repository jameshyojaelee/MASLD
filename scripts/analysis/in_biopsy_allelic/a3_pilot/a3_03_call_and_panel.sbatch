#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --array=1-22
# A3 pilot, one autosome per task:
#  (1) joint RNA genotype calls on the 30 Geuvadis BAMs over merged GENCODE v49
#      exons (the same caller settings as PRESPEC_A section 6): SNVs, MAPQ 255,
#      BQ 20, DP < 10 -> missing, QUAL >= 30, REDIportal editing sites removed;
#  (2) DNA truth: the same 30 individuals from the 1kGP high-coverage panel;
#  (3) imputation reference: the panel minus every sample sharing a FamilyID with
#      any of the 30 (biallelic SNVs only).
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out a3_03_call_and_panel.sbatch <exec_dir>
set -euo pipefail
exec_dir=$1
c=chr$SLURM_ARRAY_TASK_ID
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external
refdir=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49
panel=$ext/allelic_refs/kgp_highcov/1kGP_high_coverage_Illumina.$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
ped=$ext/allelic_refs/kgp_highcov/20130606_g1k_3202_samples_ped_population.txt
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
w=$exec_dir/work; mkdir -p "$w/calls" "$w/truth" "$w/panel" "$w/regions"

# exon regions for this chromosome
awk -F'\t' -v c="$c" '$1 == c && $3 == "exon" {print $1 "\t" $4 - 1 "\t" $5}' "$refdir/gencode.v49.primary_chr.gtf" \
  | sort -k2,2n | awk 'BEGIN {OFS = "\t"} {if (NR > 1 && $2 <= e) {if ($3 > e) e = $3} else {if (NR > 1) print c, s, e; c = $1; s = $2; e = $3}} END {print c, s, e}' \
  > "$w/regions/$c.exons.bed"
# REDIportal sites on this chromosome (header line skipped; columns 1-2 = chrom, position)
zcat "$ext/allelic_refs/rediportal/TABLE1_hg38_v3.txt.gz" | awk -F'\t' -v c="$c" 'NR > 1 && $2 == c {print $2 "\t" $3}' \
  | sort -k2,2n -u > "$w/regions/$c.rediportal.tsv"

# (1) RNA calls
awk -F'\t' 'NR > 1 {print $1}' "$ext/geuvadis/selected_30_eur.tsv" \
  | awk -v d="$ext/geuvadis/bam" '{print d "/" $1 "/" $1 ".Aligned.sortedByCoord.out.bam"}' > "$w/calls/bams.txt"
bcftools mpileup -f "$refdir/genome.primary_chr.fa" -R "$w/regions/$c.exons.bed" -q 255 -Q 20 -d 1000 \
    -a FORMAT/AD,FORMAT/DP -b "$w/calls/bams.txt" -Ou \
  | bcftools call -m -v -Ou \
  | bcftools view -v snps -m2 -M2 -i 'QUAL >= 30' -T "^$w/regions/$c.rediportal.tsv" -Ou \
  | bcftools +setGT - -Ou -- -t q -n . -i 'FMT/DP < 10' \
  | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$w/calls/$c.rna.vcf.gz"
bcftools index -t "$w/calls/$c.rna.vcf.gz"

# (2) DNA truth for the 30 individuals
awk -F'\t' 'NR > 1 {print $2}' "$ext/geuvadis/selected_30_eur.tsv" > "$w/truth/individuals.txt"
bcftools view -S "$w/truth/individuals.txt" -v snps -m2 -M2 -Oz -o "$w/truth/$c.dna.vcf.gz" "$panel"
bcftools index -t "$w/truth/$c.dna.vcf.gz"

# (3) imputation panel without the 30 individuals' families
awk 'NR == FNR {sel[$1] = 1; next} FNR == 1 {next} {row[FNR] = $0; f[FNR] = $1; s[FNR] = $2; if ($2 in sel) fam[$1] = 1}
     END {for (i in row) if (f[i] in fam) print s[i]}' "$w/truth/individuals.txt" "$ped" | sort -u > "$w/panel/excluded_samples.txt"
bcftools view -S "^$w/panel/excluded_samples.txt" -v snps -m2 -M2 -Oz -o "$w/panel/$c.ref.vcf.gz" "$panel"
bcftools index -t "$w/panel/$c.ref.vcf.gz"
echo "$c: rna sites $(bcftools view -H "$w/calls/$c.rna.vcf.gz" | wc -l); excluded panel samples $(wc -l < "$w/panel/excluded_samples.txt")"
