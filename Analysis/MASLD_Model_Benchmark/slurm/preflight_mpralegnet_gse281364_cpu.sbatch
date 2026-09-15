#!/usr/bin/env bash
#SBATCH --job-name=bench-mpralegnet-cpu-check
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
ADMISSION="${ROOT}/executions/mpralegnet-admission-21066221"
RUNTIME_PROBE="${ROOT}/executions/mpralegnet-native-probe-21066276"
COMMON_FIXTURE="${ROOT}/executions/dna-lm-sequence-fixture-21064549"
FIXTURE="${ROOT}/executions/gse281364-mpralegnet-fixture-21069059"
TARGET="${ROOT}/executions/mpralegnet-gse281364-cpu-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite MPRALegNet CPU preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af ${BASE}" \
  "d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3 ${TORCH}" \
  "5fed22c975200495ea53b38b09f5286cf13a9230a138cb72a0c758c85d5cb9d4 ${CONVERSION}" \
  "6f77ec23b4990efe416a3f69fcacb9ccce942b2fb062f17855fb9bd32e26ca7f ${ADMISSION}" \
  "afd9bfe3596e826e6987c61d3802af4ebaf2c1fe1bd67a9ab97cb62ab36c1263 ${RUNTIME_PROBE}" \
  "c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2 ${COMMON_FIXTURE}" \
  "33c3ea9fe134157b9f8f29396cd9d0d565243a9f6e44db3f458cd970c486d32c ${FIXTURE}"; do
  read -r expected directory <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${BASE}" "${TORCH}" "${CONVERSION}" "${ADMISSION}" \
  "${RUNTIME_PROBE}" "${COMMON_FIXTURE}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
install -d -m 0750 "${STAGE}/validation"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${RUNTIME_PROBE}/offline-guard:${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}"
python -m unittest -v \
  tests.unit.test_build_gse281364_mpralegnet_fixture \
  tests.unit.test_mpralegnet_predict_gse281364 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python - "${ADMISSION}" "${RUNTIME_PROBE}" "${COMMON_FIXTURE}" "${FIXTURE}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
import torch
from scripts.mpralegnet_predict_gse281364 import cpu_gpu_parity, load_model

admission, probe, common, fixture, output = map(Path, sys.argv[1:])
fixture_receipt = json.loads((fixture / "fixture/receipt.json").read_text())
if (
    fixture_receipt["status"] != "pass_outcome_blind_mpralegnet_fixture"
    or fixture_receipt["outcomes_read"]
    or fixture_receipt["reporter_counts_read"]
):
    raise SystemExit("MPRALegNet fixture preflight receipt differs")
torch.manual_seed(20260824)
torch.use_deterministic_algorithms(True)
torch.set_num_threads(4)
model = load_model(
    admission / "sources/author_model.py",
    admission / "metadata/config.json",
    admission / "checkpoint/hepg2_best_model_test1_val2.safetensors",
)
maximum = cpu_gpu_parity(
    model,
    common / "fixture/sequence_manifest.tsv",
    common / "fixture/common_6000.alleles.fa.gz",
    probe / "predictions/native_scores.npz",
)
receipt = {
    "schema_version": "masld-bench-mpralegnet-gse281364-cpu-preflight-v1",
    "status": "pass_cpu_runtime_assessment",
    "device": "cpu",
    "cpu_to_l40s_parity_max_abs_diff": maximum,
    "cpu_to_l40s_parity_tolerance": 1.0e-4,
    "cpu_prediction_eligible": maximum <= 1.0e-4,
    "required_prediction_device": "cpu" if maximum <= 1.0e-4 else "l40s",
    "checkpoint_sha256": "470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39",
    "fixture_elements": fixture_receipt["elements"],
    "outer_locus_sequence_groups": fixture_receipt["outer_locus_sequence_groups"],
    "input_lane": fixture_receipt["input_lane"],
    "outcomes_read": False,
    "reporter_counts_read": False,
    "sealed_outcomes_read": False,
    "standalone_champion_eligible": False,
}
(output / "validation/preflight_receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(receipt, sort_keys=True))
PY
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_gse281364_mpralegnet_fixture.py" \
  "${ROOT}/scripts/mpralegnet_predict_gse281364.py" \
  "${ROOT}/scripts/mpralegnet_probe_native.py" \
  "${ROOT}/tests/unit/test_build_gse281364_mpralegnet_fixture.py" \
  "${ROOT}/tests/unit/test_mpralegnet_predict_gse281364.py" \
  "${ROOT}/slurm/preflight_mpralegnet_gse281364_cpu.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/preflight_receipt.json").read_text())
if (
    receipt["status"] != "pass_cpu_runtime_assessment"
    or receipt["device"] != "cpu"
    or receipt["fixture_elements"] != 4359
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["standalone_champion_eligible"]
):
    raise SystemExit("MPRALegNet CPU preflight receipt differs")
freeze_tree(root, {
    "artifact_class": "mpralegnet_gse281364_cpu_preflight",
    "model_id": "mpralegnet_hepg2_test1_val2",
    "dataset_id": "gse281364",
    "device": "cpu",
    "cpu_prediction_eligible": receipt["cpu_prediction_eligible"],
    "required_prediction_device": receipt["required_prediction_device"],
    "outcomes_read": False,
    "standalone_champion_eligible": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'MPRALegNet GSE281364 CPU preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
