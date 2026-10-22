#!/bin/bash
#SBATCH --job-name=cgb8v2
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/logs/slurm-%x-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
TS=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/eqtl-$TS
$ENV/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/scripts/cgb8_broadaway_counts.py --out /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/eqtl-$TS/broadaway_gene_counts.tsv
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/out-$TS
$ENV/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/scripts/cgb8_gene_type.py --genes /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv --eqtl /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/eqtl-$TS/broadaway_gene_counts.tsv --out $OUT
$ENV/python -m pip freeze > $OUT/pip_freeze.txt 2>/dev/null || true
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/scripts/*.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/PRESPEC_*.md /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/DEVIATIONS.md > $OUT/script_sha256.txt
