#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00
#SBATCH --job-name=masld_nmfk --output=RNA-seq/logs/masld_nmfk_%j.out --error=RNA-seq/logs/masld_nmfk_%j.err
# Submit with: sbatch --export=ALL,NMF_K=5 RNA-seq/run_94a_single_k.sbatch
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "NMF_K=${NMF_K}"
Rscript RNA-seq/94a_nmf_single_k.R
