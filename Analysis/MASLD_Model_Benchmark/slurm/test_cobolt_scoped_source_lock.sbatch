#!/usr/bin/env bash
#SBATCH --job-name=masld-cobolt-source-scope
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
  tests.integration.test_campaign.CampaignIntegrationTests.test_scoped_cobolt_lock_ignores_queue_but_detects_bound_source \
  tests.integration.test_campaign.CampaignIntegrationTests.test_supplementary_toml_changes_registry_run_and_plan_identity \
  tests.integration.test_campaign.CampaignIntegrationTests.test_supplementary_json_changes_registry_run_and_plan_identity \
  tests.integration.test_campaign.CampaignIntegrationTests.test_source_change_changes_run_identity_but_outputs_do_not \
  tests.integration.test_campaign.CampaignIntegrationTests.test_configuration_snapshot_rejects_symlinks
