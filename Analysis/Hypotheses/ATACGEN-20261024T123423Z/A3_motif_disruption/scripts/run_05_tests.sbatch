#!/bin/bash
#SBATCH --job-name=A3_tests
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption/logs/%x_%j.err
set -euo pipefail
A3=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
$PY -s $A3/scripts/05_tests.py
