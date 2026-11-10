#!/bin/bash
#SBATCH --job-name=cgb8_posthoc
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z/logs/%x_%j.out
set -euo pipefail
H=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB8-20261022T122538Z
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$H/posthoc-$(date -u +%Y%m%dT%H%M%SZ)
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY $H/scripts/cgb8_posthoc.py --genes $R/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv \
  --eqtl $H/eqtl-20261022T125150Z/broadaway_gene_counts.tsv \
  --coloc $R/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv --out $OUT
$PY -m pip freeze > $OUT/pip_freeze.txt
sha256sum $H/scripts/*.py $H/DEVIATIONS.md > $OUT/script_sha256.txt
