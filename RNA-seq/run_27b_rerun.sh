#!/bin/bash
#SBATCH --job-name=liver_atlas_27b_rerun
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/atlas_27b_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/atlas_27b_rerun_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/

echo "=== Re-running 27b: Benchmark Presets (fixed drug target definitions) ==="
echo "Started at $(date)"
Rscript RNA-seq/27b_benchmark_presets.R
echo "=== Done at $(date) ==="
