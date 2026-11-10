#!/bin/bash
#SBATCH --job-name=gen-signal-inventory
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/inventory-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/inventory-%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
work=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics
cd "$work"
sha256sum -c SPEC.sha256
sha256sum -c inventory_code.sha256
python3 "$work/inventory_and_match.py" "$work/inventory_${SLURM_JOB_ID}"
