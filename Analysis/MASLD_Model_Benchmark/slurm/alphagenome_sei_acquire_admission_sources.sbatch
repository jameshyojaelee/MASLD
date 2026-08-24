#!/usr/bin/env bash
#SBATCH --job-name=masld-alpha-sei-admit
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/alphagenome_sei_admission.json"
TARGET="${ROOT}/executions/alphagenome-sei-admission-sources-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite AlphaGenome/Sei admission artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/sources" "${STAGE}/validation"
curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30)

ALPHA_REV=1e55dcffb98ba26b31e74edc5e9f038f54c0e89d
ALPHA_RAW="https://raw.githubusercontent.com/google-deepmind/alphagenome_research/${ALPHA_REV}"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_readme" "${ALPHA_RAW}/README.md"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_license" "${ALPHA_RAW}/LICENSE"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_pyproject" "${ALPHA_RAW}/pyproject.toml"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_dna_model" "${ALPHA_RAW}/src/alphagenome_research/model/dna_model.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_heads" "${ALPHA_RAW}/src/alphagenome_research/model/heads.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_metadata_py" "${ALPHA_RAW}/src/alphagenome_research/model/metadata/metadata.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_variant_scoring" "${ALPHA_RAW}/src/alphagenome_research/model/variant_scoring/variant_scoring.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_gene_mask" "${ALPHA_RAW}/src/alphagenome_research/model/variant_scoring/gene_mask.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_human_metadata" "${ALPHA_RAW}/src/alphagenome_research/model/metadata/OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_schemas" "${ALPHA_RAW}/src/alphagenome_research/model/schemas.py"

ALPHA_HF_REV=a8f293a76ee73d5b57f3bf2ae146510589fcf187
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_hf_model_api.json" \
  "https://huggingface.co/api/models/google/alphagenome-all-folds/revision/${ALPHA_HF_REV}"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_hf_tree_api.json" \
  "https://huggingface.co/api/models/google/alphagenome-all-folds/tree/${ALPHA_HF_REV}?recursive=true&expand=true"
curl "${curl_flags[@]}" --output "${STAGE}/sources/alphagenome_model_card.md" \
  "https://huggingface.co/google/alphagenome-all-folds/resolve/${ALPHA_HF_REV}/README.md"

SEI_REV=1e16bff9b1618af7065a23f82c4a3bb82aaeb689
SEI_RAW="https://raw.githubusercontent.com/FunctionLab/sei-framework/${SEI_REV}"
curl "${curl_flags[@]}" --output "${STAGE}/sources/sei_readme" "${SEI_RAW}/README.md"
curl "${curl_flags[@]}" --output "${STAGE}/sources/sei_license" "${SEI_RAW}/LICENSE.txt"
curl "${curl_flags[@]}" --output "${STAGE}/sources/sei_model_py" "${SEI_RAW}/model/sei.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/sei_train_yml" "${SEI_RAW}/train/train.yml"
curl "${curl_flags[@]}" --output "${STAGE}/sources/sei_variant_yml" "${SEI_RAW}/model/sei_varianteffect_prediction.yml"
curl "${curl_flags[@]}" --output "${STAGE}/sources/sei_zenodo_record_api.json" \
  "https://zenodo.org/api/records/4906997"

module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/alphagenome_sei_validate_sources.py" \
  --config "${CONFIG}" \
  --source-dir "${STAGE}/sources" \
  --project-root "${ROOT}" \
  --output "${STAGE}/review" \
  > "${STAGE}/validation/source-validation.stdout.txt" \
  2> "${STAGE}/validation/source-validation.stderr.txt"
curl --version > "${STAGE}/validation/curl.version.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/alphagenome_sei_validate_sources.py" \
  "${ROOT}/slurm/alphagenome_sei_acquire_admission_sources.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "review" / "source_admission_receipt.json").read_text())
models = receipt["model_dispositions"]
if (
    receipt["status"] != "pass"
    or receipt["checkpoint_bytes_downloaded"]
    or receipt["checkpoint_bytes_deserialized"]
    or receipt["model_forward_executed"]
    or receipt["observed_outcomes_loaded"]
    or models["alphagenome_all_folds"]["next_model_job_allowed"]
    or not models["sei"]["safe_archive_acquisition_allowed"]
    or models["sei"]["checkpoint_deserialization_allowed"]
):
    raise SystemExit("AlphaGenome/Sei source disposition differs")
freeze_tree(
    root,
    {
        "artifact_class": "alphagenome_sei_admission_sources",
        "models": ["alphagenome_all_folds", "sei"],
        "checkpoint_bytes_downloaded": False,
        "model_forward_executed": False,
        "independent_source_license_exposure_review": True,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'AlphaGenome/Sei source admission artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
