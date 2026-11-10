#!/bin/bash
#SBATCH --job-name=cga7
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGA7-20261022T191715Z/logs/slurm-%x-%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGA7-20261022T191715Z/out-$(date -u +%Y%m%dT%H%M%SZ)
$ENV/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGA7-20261022T191715Z/scripts/cga7_catalog_fields.py --out $OUT
$ENV/python -m pip freeze > $OUT/pip_freeze.txt 2>/dev/null || true
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGA7-20261022T191715Z/scripts/*.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGA7-20261022T191715Z/PRESPEC_CGA7.md > $OUT/script_sha256.txt
