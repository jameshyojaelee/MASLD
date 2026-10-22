#!/bin/bash
#SBATCH --job-name=7521
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce/jobs/7521/%j.out
set -euo pipefail
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce
FM=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $D/scripts/fig2b_render_repair.R
# Which studies have any fine-mapping row (study is column 7 canonical, 6 targeted)
awk -F, 'NR>1{c[$7]++} END{for(k in c) print "canonical\t"k"\t"c[k]}' $FM/combined_finemapping.csv > $D/out/finemapping_rows_by_study.tsv
awk -F, 'NR>1{c[$6]++} END{for(k in c) print "coloc_targeted\t"k"\t"c[k]}' $FM/combined_finemapping_coloc_targeted.csv >> $D/out/finemapping_rows_by_study.tsv
grep -c . $D/out/finemapping_rows_by_study.tsv
