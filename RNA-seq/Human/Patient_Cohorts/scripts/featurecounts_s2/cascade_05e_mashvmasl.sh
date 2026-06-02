#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_05e_mashvmasl_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_05e_mashvmasl_%j.err

# Single clean dream run (NO concurrent BiocParallel — avoids the stage1 deadlock).
# Refreshes mash_vs_masl_dream_strict.csv to -s 2 (the one contrast volcano still needs).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts
export MASH_DEF=strict

echo "=== 05e mash_vs_masl_dream_strict (single dream, -s 2) ==="
echo "Started: $(date)"
Rscript analysis/integration/scripts/05e_mash_vs_masl_dream_strict.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05e"; exit 1; fi
echo "Finished: $(date)"
