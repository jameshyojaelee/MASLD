#!/bin/bash
#SBATCH --job-name=unique-verify
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/NUCLEUS-FIX-20261107/verify_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task=Analysis/Hypotheses/NUCLEUS-FIX-20261107
rs=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
"$rs" "$task/scripts/refit_composition.R"
"$rs" "$task/scripts/verify_permutations.R"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python "$task/scripts/verify_results.py"
