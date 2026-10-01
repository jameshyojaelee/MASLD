#!/bin/bash
#SBATCH --job-name=masld_ag_resource_scope
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_resource_scope_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_resource_scope_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python3 scripts/manuscript/validate_resource_scope.py
