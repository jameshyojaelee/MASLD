#!/bin/bash
#SBATCH --job-name=V6a_test5_dmnull
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/method_comparison/test5_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/method_comparison/test5_%j.err

source ~/.bashrc
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
N_DRAWS=${N_DRAWS:-30} python -u docs/manuscript/verification/method_comparison/test5_H2H3_degree_matched.py
