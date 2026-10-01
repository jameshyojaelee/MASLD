#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex
export CODEX_REC_OUTPUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m unittest -v test_likelihood_v1 test_v12 test_v13 test_v17
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python check_predictor_workspace.py
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python check_score_path.py
