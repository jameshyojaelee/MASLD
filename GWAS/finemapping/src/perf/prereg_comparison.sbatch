#!/bin/bash
#SBATCH --job-name=coloccmp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/prereg_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/prereg_%j.err
# Pre-registered rerun-vs-canonical comparison. READ-ONLY: reads
# results/susie_coloc/ and results/susie_coloc_rerun/ and writes nothing but its log.
set -eo pipefail
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
python3 "${FM}/src/perf/prereg_comparison.py"
