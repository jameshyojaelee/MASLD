#!/bin/bash
#SBATCH --job-name=deconv
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/deconv_treat_%j.out
#SBATCH --error=logs/deconv_treat_%j.err
# ---------------------------------------------------------------------------
# Audit P2#17 remediation (2026-06-29): re-run cell-type intrinsic attribution
# with the canonical Tier-1 gate switched from raw padj<0.05 & |logFC|>0.5 (1853)
# to TREAT treat_fdr<0.05 (the 1,918 set). Output -> results/celltype_attribution/c2_treat/
# (does NOT overwrite the prior c2_recount attribution).
# ---------------------------------------------------------------------------
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
Rscript 80_celltype_intrinsic_attribution_C2.R
