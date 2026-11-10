#!/bin/bash
#SBATCH --job-name=DIG1gtest
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0-02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_rna_genotypes_20261024T235806Z/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_rna_genotypes_20261024T235806Z/%j.err
set -euo pipefail
BAM=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/GSE130970/alignments/star/SRR9036307/SRR9036307.Aligned.sortedByCoord.out.bam
BCF=/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools
FA=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa
$BCF mpileup -f $FA -R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/scripts/62_targets.tsv -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -B -d 100000 "$BAM##idx##$BAM.bai" | $BCF call -m -Oz -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_rna_genotypes_20261024T235806Z/test.vcf.gz
$BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%AD\t%DP]\n' /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/test_rna_genotypes_20261024T235806Z/test.vcf.gz
