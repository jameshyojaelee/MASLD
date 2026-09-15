#!/usr/bin/env bash
#SBATCH --job-name=build-enformer-fixture
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SPLIT="${ROOT}/executions/gse281364-outcome-blind-splits-21068885"
SEI="${ROOT}/executions/gse281364-sei-fixture-21068925"
FASTA="${ROOT}/executions/sequence-fasta-21063808"
TARGET="${ROOT}/executions/gse281364-enformer-fixture-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Enformer fixture artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "26a20a18c7697c4e16eb645640428d3fd6290957efcef4df18219f9e16320de8 ${SPLIT}/ARTIFACTS.json" \
  "620464c0ef7cf378eb2b44a1608cb6b9c273bd70b8770a43b2a95d21e725af4d ${SEI}/ARTIFACTS.json" \
  "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 ${FASTA}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}:${ROOT}/src" python -m unittest -v tests.unit.test_build_gse281364_enformer_fixture \
  > "${STAGE}/validation/tests.stdout.txt" 2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}:${ROOT}/src" python scripts/build_gse281364_enformer_fixture.py \
  --split-root "${SPLIT}" \
  --split-sha256 26a20a18c7697c4e16eb645640428d3fd6290957efcef4df18219f9e16320de8 \
  --sei-fixture-root "${SEI}" \
  --sei-fixture-sha256 620464c0ef7cf378eb2b44a1608cb6b9c273bd70b8770a43b2a95d21e725af4d \
  --fasta-root "${FASTA}" \
  --fasta-sha256 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  --output "${STAGE}/fixture" > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum scripts/build_gse281364_enformer_fixture.py \
  tests/unit/test_build_gse281364_enformer_fixture.py \
  slurm/build_gse281364_enformer_fixture.sbatch > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "fixture/receipt.json").read_text())
if (
    receipt["status"] != "pass_outcome_blind_enformer_fixture"
    or receipt["elements"] != 1033
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["outer_folds"] != 5
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["champion_eligible"]
):
    raise SystemExit("Enformer fixture receipt differs")
freeze_tree(root, {
    "artifact_class": "gse281364_outcome_blind_enformer_fixture",
    "dataset_id": "gse281364",
    "model_id": "enformer_crested_restricted_port",
    "elements": 1033,
    "outer_locus_sequence_groups": 1033,
    "outcomes_read": False,
    "champion_eligible": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 Enformer fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
