#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=2 --mem=32G --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=RNA-seq/logs/nmf_k6refit_merge_%j.out
#SBATCH --error=RNA-seq/logs/nmf_k6refit_merge_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)" && micromamba activate rnaseq && cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean_refit.rds"
export NMF_SHARD_TAG="refit"
Rscript RNA-seq/94b_merge_shards.R
