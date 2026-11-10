#!/bin/bash
#SBATCH --job-name=code-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu,io
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/pool_%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $SS/scripts/pool_sensitivity.R $SS/out_pool-$(date -u +%Y%m%dT%H%M%SZ)
echo EXIT_OK
