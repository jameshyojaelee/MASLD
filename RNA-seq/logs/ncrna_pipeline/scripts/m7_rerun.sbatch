#!/bin/bash
#SBATCH --job-name=ncRNA_M7_rerun
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "Module 7 RERUN | Start: $(date)"
Rscript scripts/figures/figS_ncrna.R
echo "Module 7 RERUN: Exit code: $? | End: $(date)"
