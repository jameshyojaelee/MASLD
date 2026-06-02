#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_step1_matrix_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_step1_matrix_%j.err
# Wave1 [A] step 1: 95_nmf_clean_ksweep.R builds the NMF input matrix from
# canonical -s2 merged_dge.rds (cache was cleared to .bak, forcing rebuild).
# Fits k=3,4 inline then quits cleanly (k=5..8 via 94a shards).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"
export NMF_RUNS=50
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== NMF step1: 95_nmf_clean_ksweep (matrix build) | $(date) ==="
Rscript RNA-seq/95_nmf_clean_ksweep.R
if [ $? -ne 0 ]; then echo "FAILED: 95_nmf_clean_ksweep"; exit 1; fi
echo "=== done $(date) ==="
