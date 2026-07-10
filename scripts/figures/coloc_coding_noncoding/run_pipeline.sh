#!/bin/bash
#SBATCH --job-name=coloc_coding_nc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --qos=interactive
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/coloc_coding_noncoding/logs/pipeline_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/coloc_coding_noncoding/logs/pipeline_%j.err
# ----------------------------------------------------------------------------
# Coding vs non-coding classification of COLOC variants (Q1 + Q2).
# Sequential 4-step pipeline. Wall time < 5 min.
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=$BASE/scripts/figures/coloc_coding_noncoding
mkdir -p $SCRIPTS/logs
cd $BASE

echo "[$(date)] [01] gather variants (A + B + C definitions)"
Rscript $SCRIPTS/01_gather_variants.R

echo "[$(date)] [02] annotate variants (VariantAnnotation + TxDb hg19)"
Rscript $SCRIPTS/02_annotate_variants.R

echo "[$(date)] [03] aggregate + DEG overlap + GWAS-ATAC overlay"
Rscript $SCRIPTS/03_aggregate_and_deg_overlap.R
echo "[$(date)] DONE"
ls -la $BASE/RNA-seq/results/coloc_variant_classes/
