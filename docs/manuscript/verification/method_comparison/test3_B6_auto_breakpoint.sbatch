#!/bin/bash
#SBATCH --job-name=V6a_test3_auto_bp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/method_comparison/test3_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/method_comparison/test3_%j.err

source ~/.bashrc
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript docs/manuscript/verification/method_comparison/test3_B6_auto_breakpoint.R
