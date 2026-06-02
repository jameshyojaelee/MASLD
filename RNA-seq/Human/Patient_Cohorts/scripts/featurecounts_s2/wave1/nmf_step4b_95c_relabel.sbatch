#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_step4b_95c_relabel_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/NMF_step4b_95c_relabel_%j.err
# Wave1 [A] step 4b: 95c relabel programs (lineage-aware) at k=6 -> program_labels.csv.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"
export NMF_CHOSEN_K=6
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== NMF step4b: 95c relabel k=6 | $(date) ==="
Rscript RNA-seq/95c_relabel_programs.R
if [ $? -ne 0 ]; then echo "FAILED: 95c_relabel_programs"; exit 1; fi
echo "=== done $(date) ==="
