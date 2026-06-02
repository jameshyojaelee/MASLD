#!/bin/bash
#SBATCH --partition=gpu --cpus-per-task=8 --mem=64G --time=48:00:00
#SBATCH --job-name=nmf_v2_k5
#SBATCH --output=RNA-seq/logs/nmf_v2_k5_%j.out
#SBATCH --error=RNA-seq/logs/nmf_v2_k5_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)" && micromamba activate rnaseq && cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export NMF_K=5
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean_v2.rds"
export NMF_SHARD_TAG="clean_v2"
Rscript RNA-seq/94a_nmf_single_k.R
