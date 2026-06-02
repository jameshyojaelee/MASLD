#!/bin/bash
#SBATCH --job-name=ggplot
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_volcano_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/s2_volcano_%j.err

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== figS_volcano_mash_masl (-s 2; all 3 contrasts now fresh) ==="
echo "Started: $(date)"
Rscript scripts/figures/figS_volcano_mash_masl.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: volcano"; exit 1; fi
echo "Finished: $(date)"
