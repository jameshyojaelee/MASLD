#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/logs/%j.out
# AN1 share-only shared-posterior check. Args: <rows.tsv> <out_dir> <chr_list|all> <mode> <n_identity>
set -euo pipefail
echo "host $(hostname) start $(date -u +%FT%TZ) args: $*"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript \
  /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/scripts/an1_share_only.R "$@"
echo "end $(date -u +%FT%TZ)"
