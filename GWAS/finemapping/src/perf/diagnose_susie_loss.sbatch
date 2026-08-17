#!/bin/bash
#SBATCH --job-name=colocdx
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/dx_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/dx_%j.err
set -eo pipefail
python3 /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/perf/diagnose_susie_loss.py
