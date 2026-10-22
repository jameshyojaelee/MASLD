#!/bin/bash
#SBATCH --job-name=6214
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/logs/%j.err
set -euo pipefail
export TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/debug/reps-20261015T130124Z/tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/scripts/02_h81_summary.R --out /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/debug/reps-20261015T130124Z --seed 20261923 --guard-dir /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z/debug/guard-20261015T130124Z/guard --debug
