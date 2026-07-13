#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=8 --mem=80G --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=RNA-seq/logs/nmf_k6refit_matrix_%j.out
#SBATCH --error=RNA-seq/logs/nmf_k6refit_matrix_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)" && micromamba activate rnaseq && cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean_refit.rds"
export NMF_SHARD_TAG="refit"
export NMF_STRIP_HLA=0
export NMF_PROTEIN_CODING_ONLY=0
export NMF_RUNS=50
Rscript RNA-seq/95_nmf_clean_ksweep.R
