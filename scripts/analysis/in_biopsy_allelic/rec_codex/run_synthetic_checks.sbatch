#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
set -euo pipefail
code=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/rec_codex
cd "$code"
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m unittest -v test_likelihood_v0
