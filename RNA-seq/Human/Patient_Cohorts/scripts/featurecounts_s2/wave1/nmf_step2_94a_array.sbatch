#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --array=3-8
#SBATCH --job-name=NMF
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_94a_k%a_%A.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_94a_k%a_%A.err
# Wave1 [A] step 2: 94a single-k NMF fit, one k per array task (k=3..8), nrun=50.
# Reads mat_nn from the clean cache built by step1; writes nmf_shard_clean_k{K}.rds.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NMF_K=$SLURM_ARRAY_TASK_ID
export NMF_RUNS=50
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"
export NMF_SHARD_TAG="clean"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== NMF 94a single-k: k=${NMF_K} nrun=${NMF_RUNS} | $(date) ==="
Rscript RNA-seq/94a_nmf_single_k.R
if [ $? -ne 0 ]; then echo "FAILED: 94a k=${NMF_K}"; exit 1; fi
echo "=== done k=${NMF_K} $(date) ==="
