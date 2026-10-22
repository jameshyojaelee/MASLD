#!/bin/bash
#SBATCH --job-name=H48_test
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H48-20261008T014825Z/dev_test_subset/test_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1 H48_NPERM=200 H48_NPERM_CAL=50
/usr/bin/time -v /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H48_mouse_knockout_ranking/h48_step2_fit.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H48-20261008T014825Z/dev_test_subset 4 test 2>&1 | grep -E "Elapsed|Maximum resident|fits done|Error|error|done" || true
/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_gears/bin/python -I /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H48_mouse_knockout_ranking/h48_step3_score.py --out /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H48-20261008T014825Z/dev_test_subset > /dev/null
ls /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H48-20261008T014825Z/dev_test_subset/step3
