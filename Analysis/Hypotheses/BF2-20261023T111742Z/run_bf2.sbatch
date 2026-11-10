#!/bin/bash
#SBATCH --job-name=bf2_order
#SBATCH --partition=io
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF2-20261023T111742Z/logs/%x_%j.out
set -euo pipefail
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF2-20261023T111742Z
OUT=$HERE/out-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p $OUT
cd $HERE && sha256sum -c PRESPEC_BF2.md.sha256
sha256sum $HERE/scripts/* > $OUT/script_sha256.txt
export OMP_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $HERE/scripts/bf2_order.R $OUT 1000
