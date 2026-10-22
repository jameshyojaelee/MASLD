#!/bin/bash
#SBATCH --job-name=resource-review
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=figures/candidates/resource-review-%j.log
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -B scripts/manuscript/build_resource_review_bundle_20261014.py
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -B scripts/manuscript/verify_resource_review_bundle_20261014.py
