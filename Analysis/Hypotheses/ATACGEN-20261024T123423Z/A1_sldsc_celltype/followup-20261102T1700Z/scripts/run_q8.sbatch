#!/bin/bash
#SBATCH --job-name=atacgen_q8
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu,io
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A1_sldsc_celltype/followup-20261102T1700Z/logs/q8_%j.log
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
F=$ROOT/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A1_sldsc_celltype/followup-20261102T1700Z
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
$RENV/bin/python -I $F/scripts/06_coloc_gene_celltype.py \
  $ROOT/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z \
  $ROOT/GWAS/finemapping/config/gwas_trait_tier.tsv \
  $ROOT/Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk \
  /gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
  $F/results/Q8
$RENV/bin/pip freeze > $F/results/Q8/pip_freeze.txt
