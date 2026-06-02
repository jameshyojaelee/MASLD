#!/bin/bash
#SBATCH --job-name=lfc_cv
#SBATCH --partition=cpu,io,bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_lfc_cv_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_lfc_cv_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== LFC threshold CV analysis (post -s 2 recount) ==="
echo "Started: $(date)"

Rscript scripts/figures/fig1g_deg_landscape_cv.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: fig1g_deg_landscape_cv.R"; exit 1; fi

echo "Finished: $(date)"
