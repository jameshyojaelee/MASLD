#!/bin/bash
#SBATCH --job-name=model-check-291
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "$ROOT"
mkdir -p executions

PYTHON="$ROOT/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
FAILED="$ROOT/executions/model-check-289-21098814.failed"
if [[ ! -e "$FAILED" ]]; then
  mkdir "$FAILED"
  cp /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/executions/model-check-289-21098814.out "$FAILED/stdout.txt"
  cp /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/executions/model-check-289-21098814.err "$FAILED/stderr.txt"
  cp "$ROOT/config/artifacts/incidents/observed_multiome_baseline_fixture_21098814.json" "$FAILED/incident.json"
fi
if [[ ! -e "$FAILED/ARTIFACTS.json" ]]; then
  PYTHONPATH="src:." "$PYTHON" -c 'from masld_bench.artifacts import freeze_tree; freeze_tree("executions/model-check-289-21098814.failed", {"artifact_class":"software_failure", "job_id":"21098814", "fit_started":False, "sealed_outcomes_accessed":False})'
fi
FAILED2="$ROOT/executions/model-check-290-21098828.failed"
if [[ ! -e "$FAILED2" ]]; then
  mkdir "$FAILED2"
  cp "$ROOT/executions/model-check-290-21098828.out" "$FAILED2/stdout.txt"
  cp "$ROOT/executions/model-check-290-21098828.err" "$FAILED2/stderr.txt"
  cp "$ROOT/config/artifacts/incidents/observed_multiome_baseline_fixture_21098828.json" "$FAILED2/incident.json"
  PYTHONPATH="src:." "$PYTHON" -c 'from masld_bench.artifacts import freeze_tree; freeze_tree("executions/model-check-290-21098828.failed", {"artifact_class":"software_failure", "job_id":"21098828", "fit_started":False, "sealed_outcomes_accessed":False})'
fi
PYTHONPATH="src:." "$PYTHON" -m unittest -v \
  tests.unit.test_observed_multiome_baselines
PYTHONPATH="src:." "$PYTHON" scripts/probe_observed_multiome_baselines.py \
  --root "$ROOT" \
  --config "$ROOT/config/campaigns/observed_multiome_baseline_fixture_20260825.json" \
  --output "$ROOT/executions/model-check-291-${SLURM_JOB_ID}"
