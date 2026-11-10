#!/bin/bash
#SBATCH --account=innovation --qos=innovation --job-name=code-review
#SBATCH --partition=cpu --prefer=v6 --cpus-per-task=8 --mem=64G --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/logs/ssx_c_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_c_nmf.R --out /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/out_c-$(date -u +%Y%m%dT%H%M%SZ)
echo EXIT_OK
