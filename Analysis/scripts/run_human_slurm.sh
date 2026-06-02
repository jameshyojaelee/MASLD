#!/bin/bash
#SBATCH --job-name=liver_human_integ
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

source /gpfs/commons/home/jameslee/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq

echo "Running Human Integration..."
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts
python3 generate_human_umap.py

echo "Done."
