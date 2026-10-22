#!/bin/bash
#SBATCH --job-name=cgb3
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB3-20261022T122538Z/logs/slurm-%x-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB3-20261022T122538Z/out-$(date -u +%Y%m%dT%H%M%SZ)
$ENV/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB3-20261022T122538Z/scripts/cgb3_genotype_triplet.py --units /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/units.tsv --out $OUT
$ENV/python -m pip freeze > $OUT/pip_freeze.txt 2>/dev/null || true
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB3-20261022T122538Z/scripts/*.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB3-20261022T122538Z/PRESPEC_*.md > $OUT/script_sha256.txt
