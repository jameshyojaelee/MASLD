#!/usr/bin/env bash
#SBATCH --job-name=model-fixture-026
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
BENCH_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT="${BENCH_ROOT}/config/gse281364_borzoi_native_screen.json"
TARGET="${BENCH_ROOT}/executions/gse281364-borzoi-native-fixture-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite native fixture artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "26a20a18c7697c4e16eb645640428d3fd6290957efcef4df18219f9e16320de8 ${BENCH_ROOT}/executions/gse281364-outcome-blind-splits-21068885/ARTIFACTS.json" \
  "620464c0ef7cf378eb2b44a1608cb6b9c273bd70b8770a43b2a95d21e725af4d ${BENCH_ROOT}/executions/gse281364-sei-fixture-21068925/ARTIFACTS.json" \
  "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 ${BENCH_ROOT}/executions/sequence-fasta-21063808/ARTIFACTS.json" \
  "c6af0a8723e5a4d8186af52c277d935c5ec4a2f4c46037b3f6954f9f77c39487 ${BENCH_ROOT}/executions/borzoi-enformer-admission-sources-21064438/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
cd "${BENCH_ROOT}"
PYTHONPATH="${BENCH_ROOT}:${BENCH_ROOT}/src" python -m unittest -v \
  tests.unit.test_build_gse281364_borzoi_native_fixture \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${BENCH_ROOT}:${BENCH_ROOT}/src" python \
  scripts/build_gse281364_borzoi_native_fixture.py \
  --project-root "${BENCH_ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/fixture" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  .gitignore \
  config/gse281364_borzoi_native_screen.json \
  scripts/build_gse281364_borzoi_native_fixture.py \
  tests/unit/test_build_gse281364_borzoi_native_fixture.py \
  slurm/build_gse281364_borzoi_native_fixture_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${BENCH_ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "fixture/receipt.json").read_text())
if (
    receipt["status"] != "pass_fixture_checkpoint_execution_blocked"
    or receipt["fixture_elements"] != 1033
    or receipt["source_locus_groups"] != 1033
    or receipt["borzoi_long_range_groups"] != 239
    or receipt["outer_folds"] != 5
    or receipt["receptive_field_buffer_bp"] != 524288
    or receipt["native_human_tracks"] != 7611
    or receipt["primary_native_tracks"] != 10
    or receipt["official_ensemble_members"] != 4
    or receipt["orientation_predictions_per_allele"] != 8
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["checkpoint_bytes_loaded"]
    or receipt["model_forward_executed"]
    or receipt["model_metrics_calculated"]
    or receipt["champion_eligible"]
):
    raise SystemExit("native fixture receipt differs")
freeze_tree(root, {
    "artifact_class": "gse281364_outcome_blind_borzoi_native_fixture",
    "dataset_id": "gse281364",
    "model_id": "borzoi_ensemble",
    "fixture_elements": 1033,
    "source_locus_groups": 1033,
    "borzoi_long_range_groups": 239,
    "receptive_field_buffer_bp": 524288,
    "native_human_tracks": 7611,
    "primary_native_tracks": 10,
    "outcomes_read": False,
    "checkpoint_bytes_loaded": False,
    "model_forward_executed": False,
    "champion_eligible": False,
    "status": "passed_fixture_checkpoint_execution_blocked",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Native fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
