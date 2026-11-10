#!/bin/bash
#SBATCH --job-name=bf4_bulk
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF4-20261024T010648Z/logs/%x_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
E=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF4-20261024T010648Z/out-20261024T011020Z
mkdir -p $O
$E/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF4-20261024T010648Z/scripts/bf4_markers.py $O
$E/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF4-20261024T010648Z/scripts/bf4_bulk.R $O $O
$E/python -m pip freeze > $O/pip_freeze.txt
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF4-20261024T010648Z/scripts/* > $O/script_sha256.txt
