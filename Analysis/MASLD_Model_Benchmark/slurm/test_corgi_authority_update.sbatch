#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-authority
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
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "${ROOT}"

printf '%s  %s\n' 8a9bdbf3dce4cbebe3815d541cfb50572c2c65212c3022e21c4fc66d66113356 \
  tests/contract/test_registry_files.py | sha256sum --check --strict
printf '%s  %s\n' 1df7633507cd716f85ab03ab7ab6a171e3e56c0e3e1f29f81b04f565134a089e \
  config/artifacts/models/corgi/checkpoints.json | sha256sum --check --strict
printf '%s  %s\n' 604c1d368770a4ab02472bd718c75b05c16e32b7f273c0f5839342561bdd35be \
  config/artifacts/models/corgi/release.json | sha256sum --check --strict
printf '%s  %s\n' cafda72a70c6f6b36c66db10f5753c3d4b268da1b64378853c39a19a4016bb0a \
  config/artifacts/models/corgi/exposure_audit.json | sha256sum --check --strict
printf '%s  %s\n' 936718cd13ffd621a4d77e2a340e9ede5f328208da4f72055a41f8c1b87af50f \
  config/artifacts/models/corgi_plus/checkpoints.json | sha256sum --check --strict
printf '%s  %s\n' 60c176303db6beb14d1b20519ada8d7aace8e2cdba94703a85c5e4428a864148 \
  config/artifacts/models/corgi_plus/exposure_audit.json | sha256sum --check --strict
printf '%s  %s\n' ae3bed391965898735121345f54eb04719129a2b3aa574f4a810719172d13ddf \
  config/models/regulatory_sequence.toml | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python -m unittest \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_roster_is_complete_unique_and_fail_closed \
  tests.contract.test_registry_files.RegistryFileContractTest.test_corgi_family_loaders_and_input_modes_fail_closed \
  -v
PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config
