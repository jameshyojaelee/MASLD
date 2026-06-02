#!/bin/bash
#SBATCH --job-name=liver_mouse_integ
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

source /gpfs/commons/home/jameslee/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq

echo "Running Mouse Integration..."
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Integration
python3 generate_mouse_umap.py

echo "Done."
