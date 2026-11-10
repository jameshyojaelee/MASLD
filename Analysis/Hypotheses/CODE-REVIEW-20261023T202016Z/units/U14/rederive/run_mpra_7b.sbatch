#!/bin/bash
#SBATCH --job-name=U14mpra
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U14/rederive/logs/%j.out
# U14: Figure 7B reporter-assay gains, paired long-range-block bootstrap.
# Runs: B 2,000 seed 9873 (seed named for 7C/7F); B 10,000 seed 20260824 (Methods TF seed);
# B 10,000 seed 20260825 (producer config seed; different generator, so not the same draws); alternates 101 and 202.
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U14/rederive/mpra_7b.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U14/rederive/out/mpra 2000 9873 10000 20260824 10000 20260825 10000 101 10000 202
$PY -m pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U14/rederive/out/mpra/pip_freeze.txt
