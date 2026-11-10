#!/bin/bash
#SBATCH --job-name=x1v2-scp3359
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/X1-20261028T103410Z
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
P=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
sha256sum -c PRESPEC.sha256
$P -I scripts/x1_xenium_v2.py out_v2
$P -m pip freeze > out_v2/pip_freeze.txt
