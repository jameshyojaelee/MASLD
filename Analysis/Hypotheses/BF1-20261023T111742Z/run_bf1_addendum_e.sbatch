#!/bin/bash
#SBATCH --job-name=bf1_addE
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF1-20261023T111742Z/logs/%x_%j.out
set -euo pipefail
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF1-20261023T111742Z
OUT=$HERE/posthoc-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p $OUT
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
cd $HERE && sha256sum -c ADDENDUM_BF1_E.md.sha256
sha256sum $HERE/scripts/bf1_addendum_e.py > $OUT/script_sha256.txt
$PY -m pip freeze > $OUT/pip_freeze.txt
$PY $HERE/scripts/bf1_addendum_e.py --out $OUT --draws 500
