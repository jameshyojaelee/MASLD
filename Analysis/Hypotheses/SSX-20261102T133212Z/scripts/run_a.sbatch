#!/bin/bash
#SBATCH --account=innovation --qos=innovation --job-name=code-review
#SBATCH --partition=cpu --prefer=v6 --cpus-per-task=4 --mem=16G --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/logs/ssx_a_%j.log
set -euo pipefail
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_a_ase_direction.R --out /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/out_a-$(date -u +%Y%m%dT%H%M%SZ)
echo EXIT_OK
