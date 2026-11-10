#!/bin/bash
#SBATCH --job-name=unique-prog-rerun
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/NUCLEUS-FIX-20261107/programs_rerun-20261109T2230Z/rerun_verify_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task=Analysis/Hypotheses/NUCLEUS-FIX-20261107
d=$task/programs_rerun-20261109T2230Z
rs=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
echo "== programs"; "$rs" "$d/rerun_programs.R"; cat "$d/COMPLETED"
echo "== composition"; "$rs" "$task/scripts/refit_composition.R"
echo "== permutation check"; "$rs" "$d/verify_permutations.R"; cat "$d/permutation_check.txt"
echo "== independent checks"; /gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python "$d/verify_results.py"
