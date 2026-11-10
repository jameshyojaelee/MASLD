#!/bin/bash
#SBATCH --account=nslab --qos=nslab --job-name=code-review
#SBATCH --partition=cpu --prefer=v6 --cpus-per-task=16 --mem=96G --time=3-00:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/logs/ssx_b_eisa_%j.log
# Usage: sbatch run_b_eisa.sbatch [time-test n genes, 0 = full run]
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TT=${1:-0}; TAG=$([ "$TT" = 0 ] && echo out_b || echo out_b_timetest)
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_b_eisa.R --out /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/${TAG}-$(date -u +%Y%m%dT%H%M%SZ) --time-test $TT
echo EXIT_OK
