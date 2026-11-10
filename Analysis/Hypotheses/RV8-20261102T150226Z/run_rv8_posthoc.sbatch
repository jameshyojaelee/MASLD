#!/bin/bash
#SBATCH --job-name=rv8ph
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV8-20261102T150226Z/logs/%x_%j.out
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV8-20261102T150226Z
sha256sum scripts/rv8_posthoc_cs_membership.py > out-20261102T150342Z/posthoc_script_sha256.txt
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -I scripts/rv8_posthoc_cs_membership.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV8-20261102T150226Z/out-20261102T150342Z
