#!/usr/bin/env bash
#SBATCH --job-name=masld-mpralegnet-admit
#SBATCH --partition=io
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
AUTHOR_REV=a20dbc7d0b08daae3af1aebcfb8a71d9f198ab5d
HF_REV=bba3b416a1450f1186b2ec3112b728650568a686
TARGET="${ROOT}/executions/mpralegnet-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed MPRALegNet admission: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite MPRALegNet admission: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/checkpoint" "${STAGE}/sources" \
  "${STAGE}/metadata" "${STAGE}/validation"
curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30 \
  --proto '=https' --tlsv1.2)
curl "${curl_flags[@]}" \
  --output "${STAGE}/checkpoint/hepg2_best_model_test1_val2.safetensors" \
  "https://huggingface.co/Ni-os/MPRALegNet/resolve/${HF_REV}/weights/hepg2_best_model_test1_val2.safetensors?download=true"
curl "${curl_flags[@]}" --output "${STAGE}/sources/author_model.py" \
  "https://raw.githubusercontent.com/autosome-ru/human_legnet/${AUTHOR_REV}/model.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/hf_model.py" \
  "https://huggingface.co/Ni-os/MPRALegNet/resolve/${HF_REV}/model.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/LICENSE" \
  "https://raw.githubusercontent.com/autosome-ru/human_legnet/${AUTHOR_REV}/LICENSE"
curl "${curl_flags[@]}" --output "${STAGE}/metadata/config.json" \
  "https://huggingface.co/Ni-os/MPRALegNet/resolve/${HF_REV}/config.json"
curl "${curl_flags[@]}" --output "${STAGE}/metadata/hepg2_config.json" \
  "https://huggingface.co/Ni-os/MPRALegNet/resolve/${HF_REV}/cell_type_configs/hepg2_config.json"
curl "${curl_flags[@]}" --output "${STAGE}/metadata/README.md" \
  "https://huggingface.co/Ni-os/MPRALegNet/resolve/${HF_REV}/README.md"

printf '%s  %s\n' \
  470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39 \
  "${STAGE}/checkpoint/hepg2_best_model_test1_val2.safetensors" \
  | sha256sum --check --strict
if [[ $(stat -c '%s' "${STAGE}/checkpoint/hepg2_best_model_test1_val2.safetensors") != 5335404 ]]; then
  printf 'MPRALegNet checkpoint size differs\n' >&2
  exit 1
fi
printf '%s  %s\n%s  %s\n%s  %s\n%s  %s\n%s  %s\n%s  %s\n' \
  8a5650b788662db07fda0e82566697b0aa1b81a1e6ff7a53448f138d00c40975 "${STAGE}/sources/author_model.py" \
  8a5650b788662db07fda0e82566697b0aa1b81a1e6ff7a53448f138d00c40975 "${STAGE}/sources/hf_model.py" \
  a43c0517b34d46afbb75b48bee74c0e52ae3612f9b860a81d566320ea49c0f15 "${STAGE}/sources/LICENSE" \
  19f9ebececfe7627f20e9ce1e5be28634f9e6c3bd5fe75847cc4457e0b5bdf09 "${STAGE}/metadata/config.json" \
  6617738400d4b442c737b9ff81f1bfe26531fd41a9d6b1cf4ee1c1d89bc3af28 "${STAGE}/metadata/hepg2_config.json" \
  d8d7a46d41a1a37fe4f0a5f637bf55c649310185329127d8a2204632e480be17 "${STAGE}/metadata/README.md" \
  | sha256sum --check --strict
cmp "${STAGE}/sources/author_model.py" "${STAGE}/sources/hf_model.py"

module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/safetensors_safe_inventory.py" \
  --input "${STAGE}/checkpoint/hepg2_best_model_test1_val2.safetensors" \
  --output "${STAGE}/validation/checkpoint_inventory.json" \
  --expected-sha256 470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39 \
  --expected-size 5335404 \
  > "${STAGE}/validation/inventory.stdout.txt"
python - "${STAGE}" <<'PY'
import ast
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
ast.parse((root / "sources/author_model.py").read_text(), filename="author_model.py")
config = json.loads((root / "metadata/config.json").read_text())
hepg2 = json.loads((root / "metadata/hepg2_config.json").read_text())
inventory = json.loads((root / "validation/checkpoint_inventory.json").read_text())
if inventory["tensor_data_loaded"] or inventory["checkpoint_deserialized"]:
    raise SystemExit("MPRALegNet safe inventory contract differs")
receipt = {
    "schema_version": "masld-bench-mpralegnet-admission-v1",
    "status": "pass",
    "checkpoint_sha256": inventory["sha256"],
    "checkpoint_size_bytes": inventory["size_bytes"],
    "tensor_count": inventory["tensor_count"],
    "parameter_and_buffer_elements": inventory["parameter_and_buffer_elements"],
    "config": config,
    "hepg2_config": hepg2,
    "author_and_hf_model_source_byte_identical": True,
    "checkpoint_provenance": "third_party_Hugging_Face_repackaging_author_equivalence_unresolved",
    "checkpoint_tensor_data_loaded": False,
    "model_source_executed": False,
    "model_instantiated": False,
    "model_forward_executed": False,
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
    "standalone_champion_eligible": False,
}
(root / "validation/admission_receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n"
)
PY
python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/safetensors_safe_inventory.py" \
  "${ROOT}/slurm/mpralegnet_acquire_admit.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/admission_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["checkpoint_sha256"] != "470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39"
    or receipt["checkpoint_size_bytes"] != 5_335_404
    or not receipt["author_and_hf_model_source_byte_identical"]
    or receipt["checkpoint_tensor_data_loaded"]
    or receipt["model_source_executed"]
    or receipt["model_instantiated"]
    or receipt["model_forward_executed"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["standalone_champion_eligible"]
):
    raise SystemExit("MPRALegNet admission receipt differs")
freeze_tree(root, {
    "artifact_class": "mpralegnet_checkpoint_admission",
    "model_id": "mpralegnet_hepg2_test1_val2",
    "checkpoint_sha256": receipt["checkpoint_sha256"],
    "code_terms": "MIT",
    "checkpoint_terms": "MIT_declared_by_model_card_provenance_equivalence_unresolved",
    "checkpoint_exposure": "target_label_unexposed",
    "checkpoint_tensor_data_loaded": False,
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
    "standalone_champion_eligible": False,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'MPRALegNet admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
