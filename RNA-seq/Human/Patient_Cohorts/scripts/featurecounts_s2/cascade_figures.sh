#!/bin/bash
#SBATCH --job-name=pub_figures
#SBATCH --partition=cpu,io,bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_figures_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_figures_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Publication Figures (post -s 2 recount) ==="
echo "Started: $(date)"

bash scripts/figures/run_pub_figures.sh 2>&1
if [ $? -ne 0 ]; then echo "WARNING: Some figures may have failed — check individual logs"; fi

echo "Finished: $(date)"
