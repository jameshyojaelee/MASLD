#!/usr/bin/env bash
#SBATCH --job-name=masld-epibert-epcot-admit
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/epibert-epcotv2-admission-sources-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
EPI_REV=dae61b434f885991c29e2fde3fbf7f70d0ff9b80
EPCOT_CODE_REV=f20362ab5139ac53448745cdca3480f952db5cf1
EPCOT_SPACE_REV=f93a2a944a49797a111f0294709416d38bedbcbe
EPCOT_DATA_REV=f445058917800811cb82c38b2c806674cc02b47f

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite admission source artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/downloads" "${STAGE}/metadata" "${STAGE}/inventories" "${STAGE}/validation"

curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30)
curl "${curl_flags[@]}" \
  --output "${STAGE}/downloads/epibert_models.tar.gz" \
  https://zenodo.org/api/records/14504049/files/epibert_models.tar.gz/content
curl "${curl_flags[@]}" \
  --output "${STAGE}/downloads/epibert_code.tar.gz" \
  "https://codeload.github.com/naumanjaved/EpiBERT/tar.gz/${EPI_REV}"
curl "${curl_flags[@]}" \
  --output "${STAGE}/downloads/epcotv2_publication_code.tar.gz" \
  "https://codeload.github.com/liu-bioinfo-lab/general_AI_model/tar.gz/${EPCOT_CODE_REV}"
curl "${curl_flags[@]}" \
  --output "${STAGE}/metadata/epibert_zenodo_14504049.json" \
  https://zenodo.org/api/records/14504049
curl "${curl_flags[@]}" \
  --output "${STAGE}/metadata/epcotv2_space_revision.json" \
  "https://huggingface.co/api/spaces/luosanj/EPCOTv2/revision/${EPCOT_SPACE_REV}"
curl "${curl_flags[@]}" \
  --output "${STAGE}/metadata/epcotv2_reference_revision.json" \
  "https://huggingface.co/api/datasets/luosanj/epcotv2_data/revision/${EPCOT_DATA_REV}"
curl "${curl_flags[@]}" \
  --output "${STAGE}/metadata/epcotv2_human_model.lfs_pointer" \
  "https://huggingface.co/spaces/luosanj/EPCOTv2/raw/${EPCOT_SPACE_REV}/models/human_model.pt"

for relative in README.md curriculum/lora_prompt_model.py func_gradio.py erna/input_region_600kb.bed; do
  destination="${STAGE}/metadata/epcotv2_space_${relative//\//_}"
  curl "${curl_flags[@]}" --output "${destination}" \
    "https://huggingface.co/spaces/luosanj/EPCOTv2/resolve/${EPCOT_SPACE_REV}/${relative}"
done

space_license_status=$(curl --location --silent --show-error \
  --output "${STAGE}/validation/epcotv2_space_LICENSE.response" \
  --write-out '%{http_code}' \
  "https://huggingface.co/spaces/luosanj/EPCOTv2/resolve/${EPCOT_SPACE_REV}/LICENSE")
reference_license_status=$(curl --location --silent --show-error \
  --output "${STAGE}/validation/epcotv2_reference_LICENSE.response" \
  --write-out '%{http_code}' \
  "https://huggingface.co/datasets/luosanj/epcotv2_data/resolve/${EPCOT_DATA_REV}/LICENSE")
if [[ "${space_license_status}" != 404 || "${reference_license_status}" != 404 ]]; then
  printf 'EPCOTv2 exact-source license status changed; reassess terms before execution\n' >&2
  exit 1
fi
rm -- "${STAGE}/validation/epcotv2_space_LICENSE.response" \
  "${STAGE}/validation/epcotv2_reference_LICENSE.response"
printf '{"epcotv2_reference_license_http_status":404,"epcotv2_space_license_http_status":404}\n' \
  > "${STAGE}/validation/undeclared_license_probe.json"

gzip -t "${STAGE}/downloads/epibert_models.tar.gz"
printf '%s  %s\n' 8c6b54ac40a5239af6a5226c4c44cbb4 \
  "${STAGE}/downloads/epibert_models.tar.gz" | md5sum --check --strict
if [[ $(stat -c '%s' "${STAGE}/downloads/epibert_models.tar.gz") != 3738627947 ]]; then
  printf 'EpiBERT model archive size differs\n' >&2
  exit 1
fi
printf '%s\n' \
  'version https://git-lfs.github.com/spec/v1' \
  'oid sha256:b824a80ed238e64e15c5a6eeae91329d24f04825c62d2a9dd328658822aba45d' \
  'size 468200299' \
  > "${STAGE}/validation/expected_epcotv2_human_model.lfs_pointer"
cmp "${STAGE}/metadata/epcotv2_human_model.lfs_pointer" \
  "${STAGE}/validation/expected_epcotv2_human_model.lfs_pointer"

module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${STAGE}/downloads/epibert_models.tar.gz" \
  --expected-size 3738627947 \
  --expected-md5 8c6b54ac40a5239af6a5226c4c44cbb4 \
  --output "${STAGE}/inventories/epibert_models.inventory.json" \
  > "${STAGE}/validation/epibert_models.inventory.stdout.txt"
python "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${STAGE}/downloads/epibert_code.tar.gz" \
  --output "${STAGE}/inventories/epibert_code.inventory.json" \
  > "${STAGE}/validation/epibert_code.inventory.stdout.txt"
python "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${STAGE}/downloads/epcotv2_publication_code.tar.gz" \
  --output "${STAGE}/inventories/epcotv2_publication_code.inventory.json" \
  > "${STAGE}/validation/epcotv2_publication_code.inventory.stdout.txt"

python - "${STAGE}" "${EPI_REV}" "${EPCOT_CODE_REV}" "${EPCOT_SPACE_REV}" "${EPCOT_DATA_REV}" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

root = Path(sys.argv[1])
epibert = json.loads((root / "inventories/epibert_code.inventory.json").read_text())
epcot = json.loads((root / "inventories/epcotv2_publication_code.inventory.json").read_text())
space = json.loads((root / "metadata/epcotv2_space_revision.json").read_text())
reference = json.loads((root / "metadata/epcotv2_reference_revision.json").read_text())
zenodo = json.loads((root / "metadata/epibert_zenodo_14504049.json").read_text())

def member_index(payload):
    return {row["path"]: row for row in payload["members"] if row["type"] == "regular_file"}

epi_files = member_index(epibert)
epcot_files = member_index(epcot)
epi_prefix = f"EpiBERT-{sys.argv[2]}/"
epcot_prefix = f"general_AI_model-{sys.argv[3]}/"
expected_epi = {
    "LICENSE": "3ee2000fdcfea66b39c931b9ee747c9c08cc183c8f3afec11b24d35987a31eda",
    "README.md": "dfa38dba7207a990f1e339532a19dc31dd03720dc3b39176d23bdc4312a19fda",
    "requirements.txt": "98f30f6c722a30659633b5e62ba4e55248bc906def576eaf6b0f9273b42f8821",
    "execute_pretraining.sh": "8e2bd3d3a268e7ffb6135ffb7158b71b6d90ef29e2a47f718e80d5cab428b9af",
    "execute_finetuning.sh": "c272415d86e874e4a807e7af00ea6d0eb6b8683ba5dca18152c3e7853882a341",
    "src/models/epibert_atac_pretrain.py": "af89aca11ddcad3364cf23a180bf3208327e56ca3cd6a6124cb9e4aa36e6e0df",
    "src/models/epibert_rampage_finetune.py": "9ae93fb8eb959fad03a152f25eaee34e061847697c4141164d9890b9f6d402b7",
    "data_processing/write_TF_records/sequence_splits/all_cell_types_map.tsv": "cee26aa3cc1f1ca322e922c2b9b814b77533e0cb085f252186f5e35c7e7b738a",
    "data_processing/write_TF_records/sequence_splits/hg38_train_524k.bed": "5ef306ef1d4395d5a4da92a601e2ed2c649b5205367f875bf51b84302ab2f3c6",
}
expected_epcot = {
    "LICENSE": "2aff8b165075fc1c015e2b0f9df0adb319e965ab51cf3806c9f99dadec797b7e",
    "README.md": "738e2190e60e5217a6076c5b014114025c280685946c746127de7a19ba89fd05",
    "requirements.txt": "2bb2625debe9305e5a88ff7decd54e6c57c237e858191834b4db81d06c16dc1c",
    "data/accession_number.xlsx": "022b45f98b625f4544407c7493c174dbc8a808987ff3225c8c4542cf2cbf999d",
    "data/input_region_600kb.bed": "d6991d59eb7469d014be9583e15d6d2799ef3af83154e4eea7648932193605d7",
    "src/model.py": "853ebbb51c4f68dbf0d6189218a18aec6961afdc83957973b049958efb37143a",
    "src/train.py": "b44ff8c6a0db1873ca339454e51fdbac80f70b199f6f4a3853ef59ff7aaeb40d",
}
for relative, digest in expected_epi.items():
    if epi_files.get(epi_prefix + relative, {}).get("sha256") != digest:
        raise SystemExit(f"EpiBERT exact source differs: {relative}")
for relative, digest in expected_epcot.items():
    if epcot_files.get(epcot_prefix + relative, {}).get("sha256") != digest:
        raise SystemExit(f"EPCOTv2 publication source differs: {relative}")
if space.get("sha") != sys.argv[4] or reference.get("sha") != sys.argv[5]:
    raise SystemExit("Hugging Face immutable revision differs")
direct_space_files = {
    "epcotv2_space_README.md": "3cd720ac2b2f6c3b145722d743ffad575a5e3af1110c4ff3a351f1f2f85ab403",
    "epcotv2_space_curriculum_lora_prompt_model.py": "2445cc49ce1e71eff25da2722df3e2812e0c75ead069d10ae784d95670804640",
    "epcotv2_space_func_gradio.py": "b5f53e0811e2e270e40b6506d863a32b16252828e506c8cf98dbec20d6a825dc",
    "epcotv2_space_erna_input_region_600kb.bed": "d6991d59eb7469d014be9583e15d6d2799ef3af83154e4eea7648932193605d7",
}
for name, digest in direct_space_files.items():
    if hashlib.sha256((root / "metadata" / name).read_bytes()).hexdigest() != digest:
        raise SystemExit(f"EPCOTv2 exact Space file differs: {name}")
license_id = ((zenodo.get("metadata") or {}).get("license") or {}).get("id")
if license_id not in {"cc-by-4.0", "CC-BY-4.0"}:
    raise SystemExit(f"EpiBERT Zenodo weight license differs: {license_id!r}")

model_inventory = json.loads((root / "inventories/epibert_models.inventory.json").read_text())
receipt = {
    "schema_version": "masld-bench-epibert-epcotv2-admission-sources-v1",
    "status": "pass",
    "epibert": {
        "code_revision": sys.argv[2],
        "code_license": "MIT",
        "weight_license": "CC-BY-4.0",
        "model_archive_sha256": model_inventory["archive_sha256"],
        "model_archive_md5": model_inventory["archive_md5"],
        "model_archive_size_bytes": model_inventory["archive_size_bytes"],
        "archive_extracted": False,
        "safe_member_inventory_passed": True,
        "terminal_disposition": "checkpoint_identity_frozen_runtime_and_numeric_parity_pending",
    },
    "epcotv2": {
        "publication_code_revision": sys.argv[3],
        "publication_code_license": "MIT",
        "paper_linked_space_revision": sys.argv[4],
        "reference_revision": sys.argv[5],
        "checkpoint_lfs_sha256": "b824a80ed238e64e15c5a6eeae91329d24f04825c62d2a9dd328658822aba45d",
        "checkpoint_size_bytes": 468200299,
        "space_code_license": "UNDECLARED",
        "weight_license": "UNDECLARED",
        "reference_license": "UNDECLARED",
        "checkpoint_object_downloaded": False,
        "checkpoint_deserialized": False,
        "terminal_disposition": "blocked_terms_no_execution",
    },
    "rna_conditioned_atac_models": [],
    "observed_atac_comparators": ["epibert", "epcotv2"],
}
(root / "admission_sources_receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY

sha256sum "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  "${ROOT}/slurm/epibert_epcotv2_acquire_admission_sources.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "admission_sources_receipt.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "observed_atac_model_admission_sources",
        "models": ["epibert", "epcotv2"],
        "epibert_model_archive_sha256": receipt["epibert"]["model_archive_sha256"],
        "epcotv2_checkpoint_object_downloaded": False,
        "epcotv2_terminal_disposition": "blocked_terms_no_execution",
        "rna_conditioned_atac_models": [],
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EpiBERT/EPCOTv2 admission sources: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
