#!/usr/bin/env bash
#SBATCH --job-name=build-model-fixture
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ENFORMER="${ROOT}/executions/gse281364-enformer-fixture-21070094"
FASTA="${ROOT}/executions/sequence-fasta-21063808"
TARGET="${ROOT}/executions/gse281364-caduceus-131k-fixture-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Caduceus fixture artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "0faf9d48511a76c282bfc412195deebda63503fb5685cd3fb2bfd77c72505231 ${ENFORMER}/ARTIFACTS.json" \
  "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 ${FASTA}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}:${ROOT}/src" python -m unittest -v \
  tests.unit.test_build_gse281364_caduceus_131k_fixture \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}:${ROOT}/src" python scripts/build_gse281364_caduceus_131k_fixture.py \
  --enformer-fixture-root "${ENFORMER}" \
  --enformer-fixture-sha256 0faf9d48511a76c282bfc412195deebda63503fb5685cd3fb2bfd77c72505231 \
  --fasta-root "${FASTA}" \
  --fasta-sha256 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  --output "${STAGE}/fixture" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  scripts/build_gse281364_caduceus_131k_fixture.py \
  tests/unit/test_build_gse281364_caduceus_131k_fixture.py \
  slurm/build_gse281364_caduceus_131k_fixture_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "fixture/receipt.json").read_text())
if (
    receipt["status"] != "pass_outcome_blind_caduceus_131k_fixture"
    or receipt["elements"] != 1033
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["outer_folds"] != 5
    or receipt["input_length_bp"] != 131072
    or receipt["pool_length_bp"] != 1536
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["checkpoint_bytes_loaded"]
    or receipt["model_forward_executed"]
    or receipt["downstream_head_fit"]
    or receipt["champion_eligible"]
):
    raise SystemExit("Caduceus 131k fixture receipt differs")
freeze_tree(root, {
    "artifact_class": "gse281364_outcome_blind_caduceus_131k_fixture",
    "dataset_id": "gse281364",
    "model_id": "caduceus",
    "elements": 1033,
    "outer_locus_sequence_groups": 1033,
    "input_length_bp": 131072,
    "outcomes_read": False,
    "checkpoint_bytes_loaded": False,
    "champion_eligible": False,
    "status": "passed",
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 Caduceus 131k fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
