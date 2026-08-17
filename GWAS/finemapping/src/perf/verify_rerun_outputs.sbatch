#!/bin/bash
#SBATCH --job-name=colocvfy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/verify_outputs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/verify_outputs_%j.err
# READ-ONLY structural verification of the rerun outputs. Writes nothing but its log.
set -eo pipefail
python3 /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/perf/verify_rerun_outputs.py
