#!/bin/bash
#SBATCH --job-name=qtllabels
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/78_qtl_labels_%j.log
set -euo pipefail
source ~/.bashrc
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src
python 78_qtl_labels.py
