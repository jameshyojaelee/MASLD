#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-freeze-208
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config
PYTHONPATH=src python -m masld_bench.cli campaign freeze \
  --campaign config/campaigns/v1_rna_atac_cobolt_5seed_smoke.toml \
  --config-root config \
  --output-root candidates
