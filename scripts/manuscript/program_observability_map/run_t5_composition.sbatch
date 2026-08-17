#!/bin/bash
#SBATCH --job-name=t5coupling
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/t5_composition_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/t5_composition_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
cd "$ROOT"

echo "=== unit tests on synthetic data with a known decomposition ==="
"$RSCRIPT" scripts/manuscript/program_observability_map/t5_tests.R

echo "=== T5 composition coupling ==="
"$RSCRIPT" scripts/manuscript/program_observability_map/t5_run_composition_coupling.R
