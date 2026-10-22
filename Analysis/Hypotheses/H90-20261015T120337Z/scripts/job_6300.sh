#!/bin/bash
#SBATCH --job-name=6300
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/logs/%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
py=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
D=Analysis/Hypotheses/H90-20261015T120337Z
"$py" $D/scripts/00_check_inputs.py
"$py" $D/scripts/01_download_gse289173.py
"$py" $D/scripts/02_bridge_discovery.py
