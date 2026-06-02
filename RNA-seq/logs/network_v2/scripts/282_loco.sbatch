#!/bin/bash
#SBATCH --job-name=netv2_282_loco
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "282 | Start: $(date)"
Rscript 282_cohort_loocv.R
echo "282 | Exit: $? | End: $(date)"
