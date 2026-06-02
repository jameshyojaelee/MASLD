#!/bin/bash
#SBATCH --job-name=figS_network_v4
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/network_v2/figS_network_v4.%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/network_v2/figS_network_v4.%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "figS_network_v4 | Start: $(date)"
Rscript scripts/figures/figS_network_v4.R
echo "figS_network_v4 | Exit: $? | End: $(date)"
