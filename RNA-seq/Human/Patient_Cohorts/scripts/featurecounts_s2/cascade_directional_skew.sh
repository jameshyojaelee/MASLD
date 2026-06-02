#!/bin/bash
#SBATCH --job-name=fgsea
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_directional_skew_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_directional_skew_%j.err

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Directional skew check (up vs down DEG enrichment + composition decomposition) ==="
echo "Started: $(date)"
Rscript RNA-seq/scripts/directional_skew_check.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: directional_skew_check.R"; exit 1; fi
echo "Finished: $(date)"
