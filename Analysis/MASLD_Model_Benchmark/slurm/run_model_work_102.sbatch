#!/usr/bin/env bash
#SBATCH --job-name=model-work-102
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
ADMISSION="${ROOT}/executions/mpralegnet-admission-21066221"
RUNTIME_PROBE="${ROOT}/executions/mpralegnet-native-probe-21066276"
COMMON_FIXTURE="${ROOT}/executions/dna-lm-sequence-fixture-21064549"
FIXTURE="${ROOT}/executions/gse281364-mpralegnet-fixture-21069059"
CPU_PREFLIGHT="${ROOT}/executions/mpralegnet-gse281364-cpu-preflight-21069356"
RUN_ID="${SLURM_JOB_ID:-preflight}"
TARGET="${ROOT}/executions/model-work-102-${RUN_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

declare -a bindings=(
  "30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
  "d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3 executions/environments/pytorch-2.6.0-cu124-portable-21064747"
  "5fed22c975200495ea53b38b09f5286cf13a9230a138cb72a0c758c85d5cb9d4 executions/environments/dna-lm-torch2.6-conversion-21064945"
  "6f77ec23b4990efe416a3f69fcacb9ccce942b2fb062f17855fb9bd32e26ca7f executions/mpralegnet-admission-21066221"
  "afd9bfe3596e826e6987c61d3802af4ebaf2c1fe1bd67a9ab97cb62ab36c1263 executions/mpralegnet-native-probe-21066276"
  "c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2 executions/dna-lm-sequence-fixture-21064549"
  "33c3ea9fe134157b9f8f29396cd9d0d565243a9f6e44db3f458cd970c486d32c executions/gse281364-mpralegnet-fixture-21069059"
  "d68c2c9a2c46e0c43f41b1b3d9d5f9e885c7345bbf11f0d54fe10685725f36da executions/mpralegnet-gse281364-cpu-preflight-21069356"
)
for binding in "${bindings[@]}"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
printf '%s  %s\n' \
  d22fcb8d87fcda9227d907d3a5e9646eca9852f3d6da5e5d16eef29703ed972f \
  "${ROOT}/scripts/mpralegnet_predict_gse281364.py" \
  | sha256sum --check --strict
printf '%s  %s\n' \
  82e87f6ef1e82ea8b3e60fa179aac30a0c21b928620f3fe7ab212e111898144d \
  "${ROOT}/tests/unit/test_mpralegnet_predict_gse281364.py" \
  | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROOT}" "${bindings[@]}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
for binding in sys.argv[2:]:
    _, relative = binding.split(" ", 1)
    verify_frozen_tree(root / relative)
preflight = json.loads(
    (root / "executions/mpralegnet-gse281364-cpu-preflight-21069356/validation/preflight_receipt.json").read_text()
)
if (
    preflight["status"] != "pass_cpu_runtime_assessment"
    or preflight["cpu_prediction_eligible"]
    or preflight["required_prediction_device"] != "l40s"
    or preflight["cpu_to_l40s_parity_max_abs_diff"] != 0.00024020671844482422
    or preflight["outcomes_read"]
    or preflight["reporter_counts_read"]
):
    raise SystemExit("frozen CPU compatibility assessment differs")
print("model-work-102 frozen inputs and L40S requirement verified")
PY
if [[ "${PREFLIGHT_ONLY:-0}" == 1 ]]; then
  printf 'model-work-102 CPU preflight passed; GPU execution remains disabled\n'
  exit 0
fi
module purge
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-work-102 artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${RUNTIME_PROBE}/offline-guard:${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}"
python - <<'PY' > "${STAGE}/validation/offline-guard.json"
import json
import socket

blocked = False
try:
    socket.socket(socket.AF_INET, socket.SOCK_STREAM)
except RuntimeError:
    blocked = True
if not blocked:
    raise SystemExit("network guard differs")
print(json.dumps({"python_INET_socket_blocked": blocked}, sort_keys=True))
PY
nvidia-smi --query-gpu=name,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi.csv"
if [[ "$(wc -l < "${STAGE}/validation/nvidia-smi.csv")" -ne 1 ]] \
  || ! grep -Fq 'NVIDIA L40S' "${STAGE}/validation/nvidia-smi.csv"; then
  printf 'Exactly one NVIDIA L40S is required.\n' >&2
  exit 1
fi
python -m unittest -v \
  tests.unit.test_build_gse281364_mpralegnet_fixture \
  tests.unit.test_mpralegnet_predict_gse281364 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python "${ROOT}/scripts/mpralegnet_predict_gse281364.py" \
  --model-source "${ADMISSION}/sources/author_model.py" \
  --config "${ADMISSION}/metadata/config.json" \
  --checkpoint "${ADMISSION}/checkpoint/hepg2_best_model_test1_val2.safetensors" \
  --manifest "${FIXTURE}/fixture/manifest.tsv" \
  --fasta "${FIXTURE}/fixture/mpralegnet.alleles.fa.gz" \
  --fixture-artifacts "${FIXTURE}/ARTIFACTS.json" \
  --probe-manifest "${COMMON_FIXTURE}/fixture/sequence_manifest.tsv" \
  --probe-fasta "${COMMON_FIXTURE}/fixture/common_6000.alleles.fa.gz" \
  --probe-scores "${RUNTIME_PROBE}/predictions/native_scores.npz" \
  --runtime-probe-artifacts "${RUNTIME_PROBE}/ARTIFACTS.json" \
  --threads 8 \
  --batch-size 512 \
  --device cuda \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/predict.stdout.json" \
  2> "${STAGE}/validation/predict.stderr.txt"
python - "${FIXTURE}" "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import math
import sys

fixture, root = map(Path, sys.argv[1:])
receipt = json.loads((root / "predictions/receipt.json").read_text())
expected_folds = {"0": 823, "1": 837, "2": 852, "3": 840, "4": 1007}
if (
    receipt["status"] != "pass_outcome_blind_prediction"
    or receipt["elements"] != 4359
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["fold_elements"] != expected_folds
    or receipt["input_length_bp"] != 230
    or receipt["forward_variant_index0"] != 115
    or receipt["reverse_complement_variant_index0"] != 114
    or receipt["cell_context"] != "HepG2"
    or receipt["input_lane"] != "GRCh38p14_230bp_genomic_window_transfer"
    or receipt["allele_effect_sign"] != "ALT_minus_REF"
    or receipt["device"] != "cuda"
    or "L40S" not in receipt["gpu"]
    or receipt["runtime_to_frozen_l40s_parity_max_abs_diff"] > 1.0e-4
    or receipt["deterministic_repeat_max_abs_diff"] != 0.0
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["model_fitted_or_adapted"]
    or receipt["standalone_champion_eligible"]
):
    raise SystemExit("MPRALegNet prediction receipt differs")
with (fixture / "fixture/manifest.tsv").open(newline="") as handle:
    fixture_rows = list(csv.DictReader(handle, delimiter="\t"))
with (root / "predictions/allele_scores.tsv").open(newline="") as handle:
    score_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(score_rows) != 4359 or [row["element_id"] for row in score_rows] != [row["element_id"] for row in fixture_rows]:
    raise SystemExit("prediction row identity differs")
for observed, expected in zip(score_rows, fixture_rows, strict=True):
    if (
        observed["outer_locus_sequence_group_id"] != expected["outer_locus_sequence_group_id"]
        or observed["outer_fold"] != expected["outer_fold"]
        or observed["ref"] != expected["ref"]
        or observed["alt"] != expected["alt"]
        or not math.isclose(
            float(observed["alt_minus_ref_score"]),
            float(observed["alternative_orientation_mean_score"])
            - float(observed["reference_orientation_mean_score"]),
            rel_tol=0.0,
            abs_tol=1.0e-7,
        )
    ):
        raise SystemExit("prediction orientation, grouping, or allele sign differs")
(root / "validation/frozen_validation.json").write_text(
    json.dumps(
        {
            "status": "pass_outcome_blind_l40s_prediction_validation",
            "elements": 4359,
            "outer_locus_sequence_groups": 1033,
            "fold_elements": expected_folds,
            "device": "l40s",
            "outcomes_read": False,
            "reporter_counts_read": False,
            "standalone_champion_eligible": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
PY
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_gse281364_mpralegnet_fixture.py" \
  "${ROOT}/scripts/mpralegnet_predict_gse281364.py" \
  "${ROOT}/scripts/mpralegnet_probe_native.py" \
  "${ROOT}/tests/unit/test_build_gse281364_mpralegnet_fixture.py" \
  "${ROOT}/tests/unit/test_mpralegnet_predict_gse281364.py" \
  "${ROOT}/slurm/run_model_work_102.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
validation = json.loads((root / "validation/frozen_validation.json").read_text())
if (
    validation["status"] != "pass_outcome_blind_l40s_prediction_validation"
    or validation["elements"] != 4359
    or validation["outer_locus_sequence_groups"] != 1033
    or validation["outcomes_read"]
    or validation["reporter_counts_read"]
    or validation["standalone_champion_eligible"]
):
    raise SystemExit("MPRALegNet frozen validation differs")
freeze_tree(
    root,
    {
        "artifact_class": "outcome_blind_MPRALegNet_HepG2_reporter_score_comparator",
        "dataset_id": "gse281364",
        "model_id": "mpralegnet_hepg2_test1_val2",
        "device": "l40s",
        "elements": 4359,
        "outer_locus_sequence_groups": 1033,
        "input_lane": "GRCh38p14_230bp_genomic_window_transfer",
        "assay_insert_input_status": "unsupported_107bp_not_padded_or_truncated",
        "allele_effect_sign": "ALT_minus_REF",
        "outcomes_read": False,
        "model_fitted_or_adapted": False,
        "standalone_champion_eligible": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-work-102 artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
