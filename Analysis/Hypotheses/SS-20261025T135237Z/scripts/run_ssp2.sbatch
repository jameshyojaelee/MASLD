#!/bin/bash
#SBATCH --job-name=code-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu,io
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
# SS SSP2 worker: --export=ALL,O=<out>,MODE=fit|summarize,I=<i>,K=<k>
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/ssp2_control_nas.R
if [ "$MODE" = fit ]; then $R $S fit "$O" "$I" "$K"; else $R $S summarize "$O"; fi
echo EXIT_OK
