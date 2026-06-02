#!/bin/bash
#SBATCH --job-name=atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/atlas_27ab_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/atlas_27ab_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/

echo "=== Running 27a: Assemble Evidence Atlas ==="
echo "Started at $(date)"
Rscript RNA-seq/27a_assemble_evidence_atlas.R
echo ""
echo "=== Running 27b: Benchmark Presets ==="
echo "Started at $(date)"
Rscript RNA-seq/27b_benchmark_presets.R
echo ""
echo "=== Done at $(date) ==="
