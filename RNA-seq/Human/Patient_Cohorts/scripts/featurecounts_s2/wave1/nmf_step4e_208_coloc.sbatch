#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_step4e_208_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_step4e_208_coloc_%j.err
# Wave1 [A] step 4e (ATLAS FEEDER): 208 subtype-COLOC -> program_causal_scores.csv (read by 217).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"
export NMF_CHOSEN_K=6
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== NMF step4e: 208 subtype COLOC | $(date) ==="
Rscript RNA-seq/208_subtype_coloc.R
if [ $? -ne 0 ]; then echo "FAILED: 208_subtype_coloc"; exit 1; fi
echo "=== done $(date) ==="
