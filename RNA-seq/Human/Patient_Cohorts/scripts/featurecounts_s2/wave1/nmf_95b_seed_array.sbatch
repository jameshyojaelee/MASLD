#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --array=0-1
#SBATCH --job-name=NMF
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_95b_seed_k%a_%A.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_95b_seed_k%a_%A.err
# Wave1 [A] step 5 (TERMINAL): 95b seed-sensitivity. One k per task: task0=k4, task1=k6.
# Seeds {42,100,2024} x nrun=50 (hardcoded in 95b). Reads mat_nn from clean cache.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
K_VALUES=(4 6)
export NMF_K_LIST="${K_VALUES[$SLURM_ARRAY_TASK_ID]}"
export NMF_RUNS=50
export NMF_REUSE_SEED42_CACHE=FALSE
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== NMF 95b seed sensitivity: k=${NMF_K_LIST} seeds={42,100,2024} | $(date) ==="
Rscript RNA-seq/95b_nmf_seed_sensitivity.R
if [ $? -ne 0 ]; then echo "FAILED: 95b k=${NMF_K_LIST}"; exit 1; fi
echo "=== done k=${NMF_K_LIST} $(date) ==="
