#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=fig5_convergence
#SBATCH --output=logs/fig5_convergence_%j.out
#SBATCH --error=logs/fig5_convergence_%j.err

# Rewritten 2026-04-22: 6-fig → 5-fig restructure. Legacy fig6_therapeutic_windows.R
# content absorbed into fig5_convergence.R; the fig6 script lives in
# archive/fig6_ditched_2026-04-22/.

source ~/.bashrc
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/figures/fig5_convergence.R
