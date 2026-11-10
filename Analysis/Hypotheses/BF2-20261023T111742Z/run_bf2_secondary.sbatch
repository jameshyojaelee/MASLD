#!/bin/bash
#SBATCH --job-name=bf2_sec
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF2-20261023T111742Z/logs/%x_%j.out
set -euo pipefail
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF2-20261023T111742Z
OUT=$HERE/secondary-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p $OUT
sha256sum $HERE/scripts/bf2_secondary.py > $OUT/script_sha256.txt
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY -m pip freeze > $OUT/pip_freeze.txt
$PY $HERE/scripts/bf2_secondary.py $HERE/out-20261023T112740Z $OUT
