#!/bin/bash
#SBATCH --job-name=A2_extract
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 4
#SBATCH --mem=110G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage/logs/extract_%j.out
set -euo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 TMPDIR=$A/work
E=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
$E/python $A/scripts/01_extract.py
$E/Rscript $A/scripts/02_liftover.R
$E/python -m pip freeze > $A/pip_freeze.txt
