#!/usr/bin/env bash
#SBATCH --job-name=agp_f2-haplotype-readout
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-readout-%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export F2_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python f2_07_deposit_readout_check.py
echo F2_READOUT_DONE
