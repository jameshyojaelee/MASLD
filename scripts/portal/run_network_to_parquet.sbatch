#!/bin/bash
#SBATCH --job-name=portal
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/logs/network_to_parquet_%j.log

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$ROOT"
mkdir -p "$ROOT/scripts/portal/logs"
cd "$ROOT"
micromamba run -n spatial python scripts/portal/network_to_parquet.py
