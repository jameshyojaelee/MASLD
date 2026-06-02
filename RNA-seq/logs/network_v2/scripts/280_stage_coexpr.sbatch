#!/bin/bash
#SBATCH --job-name=netv2_280_stage_coexpr
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "280 | Start: $(date)"
Rscript 280_stage_stratified_coexpr.R
echo "280 | Exit: $? | End: $(date)"
