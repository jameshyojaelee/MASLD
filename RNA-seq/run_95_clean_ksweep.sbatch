#!/bin/bash
#SBATCH --partition=bigmem --cpus-per-task=16 --mem=220G --time=48:00:00
#SBATCH --job-name=masld_nmfclean --output=RNA-seq/logs/masld_nmfclean_%j.out --error=RNA-seq/logs/masld_nmfclean_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
# Optional env vars (default off):
#   NMF_STRIP_HLA=1
#   NMF_PROTEIN_CODING_ONLY=1
Rscript RNA-seq/95_nmf_clean_ksweep.R
