#!/usr/bin/env bash
#SBATCH --job-name=masld-artifact-publication
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1

module load python/3.11.5-GCCcore-13.2.0
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

PYTHONPATH=src python -m unittest tests.unit.test_artifacts -v
candidate_root="candidates/fallback-validation-${SLURM_JOB_ID}"
freeze_json="$(PYTHONPATH=src python -m masld_bench.cli campaign freeze \
  --campaign config/campaigns/v1_admission.toml \
  --config-root config \
  --output-root "${candidate_root}")"
printf '%s\n' "${freeze_json}"
candidate_path="$(python -c 'import json,sys; print(json.load(sys.stdin)["candidate"])' <<<"${freeze_json}")"
PYTHONPATH=src python -m masld_bench.cli campaign verify \
  --candidate "${candidate_path}"
