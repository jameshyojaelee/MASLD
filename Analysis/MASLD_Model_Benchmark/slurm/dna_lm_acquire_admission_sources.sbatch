#!/usr/bin/env bash
#SBATCH --job-name=masld-dna-lm-admit
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
CONFIG="${ROOT}/config/dna_lm_admission.json"
TARGET="${ROOT}/executions/dna-lm-admission-sources-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite DNA-LM source artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/sources" "${STAGE}/validation"
curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30)

module load python/3.11.5-GCCcore-13.2.0
python - "${CONFIG}" <<'PY' |
import json
from pathlib import Path
import sys

config = json.loads(Path(sys.argv[1]).read_text())
for name, record in config["source_files"].items():
    if not name.replace("_", "").isalnum() or "\t" in record["url"]:
        raise SystemExit("source name or URL differs")
    print(f"{name}\t{record['url']}")
PY
while IFS=$'\t' read -r name url; do
  curl "${curl_flags[@]}" --output "${STAGE}/sources/${name}" "${url}"
done

python - "${CONFIG}" <<'PY' |
import json
from pathlib import Path
import sys

config = json.loads(Path(sys.argv[1]).read_text())
for model_id, model in config["models"].items():
    print(f"{model_id}\t{model['checkpoint_repository']}\t{model['checkpoint_revision']}")
PY
while IFS=$'\t' read -r model_id repository revision; do
  curl "${curl_flags[@]}" \
    --output "${STAGE}/sources/${model_id}_hf_model_api.json" \
    "https://huggingface.co/api/models/${repository}/revision/${revision}"
  curl "${curl_flags[@]}" \
    --output "${STAGE}/sources/${model_id}_hf_tree_api.json" \
    "https://huggingface.co/api/models/${repository}/tree/${revision}?recursive=true&expand=true"
done

python "${ROOT}/scripts/dna_lm_validate_sources.py" \
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
  "${ROOT}/scripts/dna_lm_validate_sources.py" \
  "${ROOT}/slurm/dna_lm_acquire_admission_sources.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "review" / "source_admission_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["checkpoint_bytes_downloaded"]
    or receipt["checkpoint_bytes_deserialized"]
    or receipt["model_forward_executed"]
    or receipt["observed_outcomes_loaded"]
    or any(
        not row["safe_checkpoint_download_allowed"]
        or row["checkpoint_deserialization_allowed"]
        or row["gpu_forward_allowed"]
        or row["embeddings_biologically_interpretable_without_head"]
        or row["supervised_eqtl_or_ieqtl_head_allowed"]
        for row in receipt["model_dispositions"].values()
    )
):
    raise SystemExit("DNA-LM source disposition differs")
freeze_tree(
    root,
    {
        "artifact_class": "dna_lm_admission_sources",
        "models": [
            "dnabert2",
            "nucleotide_transformer",
            "hyenadna",
            "caduceus",
            "evo",
        ],
        "checkpoint_bytes_downloaded": False,
        "model_forward_executed": False,
        "independent_source_license_exposure_review": True,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'DNA-LM source admission artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
