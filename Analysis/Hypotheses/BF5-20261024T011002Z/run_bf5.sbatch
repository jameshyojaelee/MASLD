#!/bin/bash
#SBATCH --job-name=bf5
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF5-20261024T011002Z/logs/%x_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
P=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF5-20261024T011002Z/out-20261024T011301Z
mkdir -p $O
$P /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF5-20261024T011002Z/scripts/bf5_spatial_control.py --out $O
$P /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF5-20261024T011002Z/scripts/bf5_nuclei.py --out $O
$P -m pip freeze > $O/pip_freeze.txt
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF5-20261024T011002Z/scripts/*.py > $O/script_sha256.txt
