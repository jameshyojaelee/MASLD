#!/bin/bash
#SBATCH --job-name=bf3_build
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z/logs/%x_%j.out
set -euo pipefail
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z/build-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p $OUT
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z && sha256sum -c PRESPEC_BF3.md.sha256
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z/scripts/bf3_build.py > $OUT/script_sha256.txt
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY -m pip freeze > $OUT/pip_freeze.txt
$PY /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z/scripts/bf3_build.py $OUT
