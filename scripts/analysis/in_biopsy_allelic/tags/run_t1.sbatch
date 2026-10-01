#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-22
# T1 per autosome: exonic tag SNVs for GTEx v8 liver eQTL leads (EUR and EAS LD).
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out run_t1.sbatch <out_dir> [exclude_related.txt]
set -euo pipefail
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/allelic_refs
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
c=chr$SLURM_ARRAY_TASK_ID
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/tags/t1_tags.py \
  --chrom "$c" --egenes "$ext/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz" \
  --gtf /gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf \
  --panel "$ext/kgp_highcov/1kGP_high_coverage_Illumina.$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" \
  --ped "$ext/kgp_highcov/20130606_g1k_3202_samples_ped_population.txt" --out "$1/tags" ${2:+--exclude-samples "$2"}
