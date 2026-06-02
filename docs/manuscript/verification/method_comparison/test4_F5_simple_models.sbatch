#!/bin/bash
#SBATCH --job-name=V6a_test4_simple
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/method_comparison/test4_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/method_comparison/test4_%j.err

source ~/.bashrc
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python docs/manuscript/verification/method_comparison/test4_F5_simple_models.py
