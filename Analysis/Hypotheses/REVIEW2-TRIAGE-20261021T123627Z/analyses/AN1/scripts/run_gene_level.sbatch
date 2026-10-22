#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/logs/%j.out
# AN1 gene-level sensitivity. Args: --out <new dir> --shares <row_shares.tsv> [...]
set -euo pipefail
echo "host $(hostname) start $(date -u +%FT%TZ) args: $*"
env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/scripts/an1_gene_level.py "$@"
echo "end $(date -u +%FT%TZ)"
