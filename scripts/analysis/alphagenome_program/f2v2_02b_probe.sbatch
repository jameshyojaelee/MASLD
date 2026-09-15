#!/usr/bin/env bash
#SBATCH --job-name=agp_f2-haplotype-v2-probe
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-v2-probe-%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export F2V2_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/f2-haplotype-v2-20260915T063000Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
"$PY" f2v2_02b_loader_equivalence.py
echo F2V2_PROBE_DONE
