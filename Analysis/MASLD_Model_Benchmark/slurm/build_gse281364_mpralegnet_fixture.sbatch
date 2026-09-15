#!/usr/bin/env bash
#SBATCH --job-name=bench-mpralegnet-fixture
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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
SPLIT_SHA=26a20a18c7697c4e16eb645640428d3fd6290957efcef4df18219f9e16320de8
REFERENCE="${ROOT}/executions/sequence-fasta-21063808"
REFERENCE_SHA=7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39
TARGET="${ROOT}/executions/gse281364-mpralegnet-fixture-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite MPRALegNet fixture artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in "${SPLIT_SHA} ${SPLIT}/ARTIFACTS.json" "${REFERENCE_SHA} ${REFERENCE}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SPLIT}" "${REFERENCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
install -d -m 0750 "${STAGE}/validation"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_build_gse281364_mpralegnet_fixture \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" python "${ROOT}/scripts/build_gse281364_mpralegnet_fixture.py" \
  --split-root "${SPLIT}" --split-sha256 "${SPLIT_SHA}" \
  --reference-root "${REFERENCE}" --reference-sha256 "${REFERENCE_SHA}" \
  --output "${STAGE}/fixture" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_gse281364_mpralegnet_fixture.py" \
  "${ROOT}/scripts/alphagenome_sei_build_fixture.py" \
  "${ROOT}/tests/unit/test_build_gse281364_mpralegnet_fixture.py" \
  "${ROOT}/slurm/build_gse281364_mpralegnet_fixture.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "fixture/receipt.json").read_text())
if (
    receipt["status"] != "pass_outcome_blind_mpralegnet_fixture"
    or receipt["elements"] != 4359
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["sequences"] != 17436
    or receipt["input_length_bp"] != 230
    or receipt["forward_variant_index0"] != 115
    or receipt["reverse_complement_variant_index0"] != 114
    or receipt["cell_context"] != "HepG2"
    or receipt["input_lane"] != "GRCh38p14_230bp_genomic_window_transfer"
    or receipt["assay_insert_input_status"] != "unsupported_107bp_not_padded_or_truncated"
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["standalone_champion_eligible"]
):
    raise SystemExit("GSE281364 MPRALegNet fixture receipt differs")
freeze_tree(root, {
    "artifact_class": "gse281364_outcome_blind_mpralegnet_fixture",
    "dataset_id": "gse281364",
    "model_id": "mpralegnet_hepg2_test1_val2",
    "elements": 4359,
    "outer_locus_sequence_groups": 1033,
    "input_lane": receipt["input_lane"],
    "outcomes_read": False,
    "standalone_champion_eligible": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 MPRALegNet fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
