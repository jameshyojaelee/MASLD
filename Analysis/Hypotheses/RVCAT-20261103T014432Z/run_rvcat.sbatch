#!/bin/bash
#SBATCH --job-name=rvcat
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RVCAT-20261103T014432Z/logs/%x_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RVCAT-20261103T014432Z/out-$(date -u +%Y%m%dT%H%M%SZ)
$PY -I /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RVCAT-20261103T014432Z/scripts/rvcat_build.py $O
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RVCAT-20261103T014432Z/scripts/rvcat_build.py > $O/script_sha256.txt
$PY -m pip freeze > $O/pip_freeze.txt
