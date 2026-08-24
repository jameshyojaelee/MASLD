#!/usr/bin/env bash
#SBATCH --job-name=masld-bench-action-scope
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
  tests.unit.test_adapters \
  tests.integration.test_campaign.CampaignIntegrationTests.test_scientific_run_requires_and_executes_frozen_actions \
  tests.integration.test_campaign.CampaignIntegrationTests.test_smoke_view_activates_only_the_frozen_parent_view \
  tests.contract.test_registry_files.RegistryFileContractTest.test_first_geneformer_smoke_lane_is_exact_and_fail_closed
