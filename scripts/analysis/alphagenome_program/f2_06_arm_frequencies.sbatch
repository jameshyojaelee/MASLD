#!/usr/bin/env bash
#SBATCH --job-name=agp_f2-haplotype-armfreq
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-armfreq-%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export F2_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python f2_06_arm_frequencies.py
echo F2_ARMFREQ_DONE
