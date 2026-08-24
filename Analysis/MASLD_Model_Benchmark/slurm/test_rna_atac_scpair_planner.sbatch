#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-scpair-plan-test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1

module load python/3.11.5-GCCcore-13.2.0
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

PYTHONPATH=src python -m unittest -v \
  tests.integration.test_campaign.CampaignIntegrationTests.test_rna_atac_reference_and_baseline_gates_use_runnable_models \
  tests.integration.test_campaign.CampaignIntegrationTests.test_model_hyperparameters_replace_instead_of_merge \
  tests.contract.test_registry_files.RegistryFileContractTest.test_every_repo_path_anchored_audit_hash_rederives \
  tests.contract.test_registry_files.RegistryFileContractTest.test_atac_profile_and_latent_capabilities_are_disjoint
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config >/dev/null
