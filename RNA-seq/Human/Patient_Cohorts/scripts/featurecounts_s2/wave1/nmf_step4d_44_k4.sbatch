#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_step4d_44_k4_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_step4d_44_k4_%j.err
# Wave1 [A] step 4d: 44 molecular subtyping at NMF_CHOSEN_K=4 (secondary k).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"
export NMF_CHOSEN_K=4
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== NMF step4d: 44 molecular subtyping k=4 | $(date) ==="
Rscript RNA-seq/44_molecular_subtyping.R
if [ $? -ne 0 ]; then echo "FAILED: 44 k=4"; exit 1; fi
echo "=== done $(date) ==="
