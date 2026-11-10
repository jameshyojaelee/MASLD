#!/bin/bash
#SBATCH --job-name=A2_one
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage/logs/one_%j.out
# Usage: sbatch run_one_test.sbatch <script stem>   debug run of one step into results_test/
set -euo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 TMPDIR=$A/work A2_ALLOW_NO_LD=1 A2_RES=$A/results_test
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python $A/scripts/$1.py
