#!/bin/bash
#SBATCH --job-name=cg0_magma
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/logs/slurm-%x-%j.out
set -euo pipefail
bash /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/scripts/cg0_magma.sh /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/magma/run-$(date -u +%Y%m%dT%H%M%SZ)
