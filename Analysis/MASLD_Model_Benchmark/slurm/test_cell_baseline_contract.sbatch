#!/usr/bin/env bash
#SBATCH --job-name=masld-cell-baseline-contract
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1

module load python/3.11.5-GCCcore-13.2.0
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

PYTHONPATH=src python -m unittest -v \
  tests.unit.test_hvg_pca_logistic \
  tests.contract.test_registry_files.RegistryFileContractTest.test_first_real_cell_baseline_is_exact_and_submission_guarded \
  tests.contract.test_registry_files.RegistryFileContractTest.test_classical_cell_baseline_screen_is_exact_and_submission_guarded \
  tests.integration.test_campaign.CampaignIntegrationTests.test_incremental_smoke_can_run_one_explicit_mandatory_baseline \
  tests.integration.test_campaign.CampaignIntegrationTests.test_incremental_smoke_rejects_nonbaseline_only_request \
  tests.integration.test_campaign.CampaignIntegrationTests.test_production_submission_is_single_use_and_records_job_id
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config
freeze_json="$(PYTHONPATH=src python -m masld_bench.cli campaign freeze \
  --campaign config/campaigns/v1_cell_classical_baselines_smoke.toml \
  --config-root config \
  --output-root "candidates/baseline-validation-${SLURM_JOB_ID}")"
printf '%s\n' "${freeze_json}"
candidate_path="$(python -c 'import json,sys; print(json.load(sys.stdin)["candidate"])' <<<"${freeze_json}")"
PYTHONPATH=src python -m masld_bench.cli campaign verify --candidate "${candidate_path}"
