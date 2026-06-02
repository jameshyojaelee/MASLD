#!/bin/bash
#SBATCH --job-name=stg117_bifurcation_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/logs/117_bifurcation_%j.out
#SBATCH --error=scripts/logs/117_bifurcation_%j.err

source ~/.bashrc
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

python scripts/117_bifurcation_divergence.py
