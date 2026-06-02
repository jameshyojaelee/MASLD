#!/bin/bash
#SBATCH --job-name=netv2_293_sankey
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "293 | Start: $(date)"
Rscript scripts/figures/figS_network_f2_switch.R
echo "293 | Exit: $? | End: $(date)"
