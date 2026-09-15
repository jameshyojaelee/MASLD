#!/usr/bin/env bash
#SBATCH --job-name=tf-tests
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
python="${root}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
target="${root}/executions/transcriptformer-script-tests-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer script tests.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}"
cd "${root}"
PYTHONPATH="${root}" "${python}" -m unittest \
  tests.unit.test_build_transcriptformer_unlabeled_fixture \
  tests.contract.test_registry_files.RegistryFileContractTest.test_transcriptformer_preflight_fails_closed_on_exposure \
  > "${stage}/unittest.txt" 2>&1
"${python}" -m py_compile \
  scripts/build_transcriptformer_unlabeled_fixture.py \
  scripts/preflight_transcriptformer_cpu.py \
  scripts/transcriptformer_extract_unlabeled_embeddings.py
PYTHONPATH="${root}/src" "${python}" - <<'PY' > "${stage}/activation-tree.json"
from pathlib import Path
import json
from masld_bench.artifacts import verify_frozen_tree

artifact = Path("executions/model-training-21071475").resolve(strict=True)
manifest = verify_frozen_tree(artifact)
metadata = manifest["metadata"]
if (
    metadata["status"] != "passed"
    or metadata["sealed_outcomes_read"]
    or metadata["sealed_champion_eligible"]
):
    raise SystemExit("TranscriptFormer activation artifact firewall differs")
print(json.dumps(metadata, sort_keys=True))
PY
PYTHONPATH="${root}/src" "${python}" - <<'PY' > "${stage}/evaluation-tree.json"
from pathlib import Path
import json
from masld_bench.artifacts import verify_frozen_tree

artifact = Path("executions/model-evaluation-21073148").resolve(strict=True)
manifest = verify_frozen_tree(artifact)
metadata = manifest["metadata"]
if (
    metadata["status"] != "passed"
    or not metadata["development_only"]
    or metadata["sealed_outcomes_read"]
    or metadata["sealed_champion_eligible"]
):
    raise SystemExit("TranscriptFormer evaluation artifact firewall differs")
print(json.dumps(metadata, sort_keys=True))
PY
"${python}" -m pip freeze --all > "${stage}/pip-freeze.txt"
sha256sum \
  scripts/build_transcriptformer_unlabeled_fixture.py \
  scripts/preflight_transcriptformer_cpu.py \
  scripts/transcriptformer_extract_unlabeled_embeddings.py \
  config/adaptation/transcriptformer.toml \
  config/artifacts/models/transcriptformer/checkpoints.json \
  config/artifacts/models/transcriptformer/development_disposition.json \
  config/artifacts/models/transcriptformer/exposure_audit.json \
  tests/unit/test_build_transcriptformer_unlabeled_fixture.py \
  tests/contract/test_registry_files.py \
  slurm/test_transcriptformer_activation_scripts.sbatch \
  > "${stage}/source.sha256"
(
  cd "${stage}"
  sha256sum activation-tree.json evaluation-tree.json pip-freeze.txt source.sha256 unittest.txt > SHA256SUMS
)
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'TranscriptFormer activation script tests: %s\n' "${target}"
