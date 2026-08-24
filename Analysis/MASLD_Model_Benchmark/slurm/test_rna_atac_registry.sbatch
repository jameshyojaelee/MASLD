#!/usr/bin/env bash
#SBATCH --job-name=masld-bench-rna-atac-registry
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
  tests.contract.test_registry_files.RegistryFileContractTest.test_every_repo_path_anchored_audit_hash_rederives \
  tests.contract.test_registry_files.RegistryFileContractTest.test_graph_and_transparent_baseline_audits_freeze_exact_boundaries \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_roster_is_complete_unique_and_fail_closed \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_census_maps_exactly_to_audit_bundles \
  tests.contract.test_registry_files.RegistryFileContractTest.test_every_campaign_binds_the_frozen_model_census \
  tests.contract.test_registry_files.RegistryFileContractTest.test_every_model_audit_bundle_is_complete_and_hash_bound \
  tests.contract.test_registry_files.RegistryFileContractTest.test_atac_profile_and_latent_capabilities_are_disjoint
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config
