#!/usr/bin/env bash
#SBATCH --job-name=model-runtime-042
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE_ENV=/gpfs/commons/home/jameslee/micromamba/envs/decima
AUDIT=${ROOT}/executions/decima-existing-cellfm-runtime-audit-21066298
BUNDLE=${ROOT}/executions/geneformer-acquisition-21014827/bundles/geneformer_v2_316m--11f311a5078fc4c8
INVENTORY=${ROOT}/executions/geneformer_v2_316m-safe-inventory-21066211/checkpoint_inventory.json
TOKENS=${ROOT}/executions/geneformer-sweep/21015001/geneformer_v2_316m__two_layer_mlp/prepare/rank_value_tokens.jsonl
TOKEN_SUMMARY=${ROOT}/executions/geneformer-sweep/21015001/geneformer_v2_316m__two_layer_mlp/prepare/tokenization_summary.json
PRIOR_EMBEDDINGS=${ROOT}/executions/embeddings/geneformer_v2_316m/embeddings.npy
PRIOR_ROW_ORDER=${ROOT}/executions/embeddings/geneformer_v2_316m/embedding_row_order.txt
ADAPTER=${ROOT}/src/masld_bench/adapters/geneformer_common_lane.py
EXTRACTOR=${ROOT}/executions/fiveseed/extract_geneformer_emb.py
TARGET=${ROOT}/executions/geneformer-v2-316m-compat-runtime-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
WORK=$(mktemp -d "${TMPDIR:-/tmp}/geneformer-runtime.${SLURM_JOB_ID}.XXXXXX")

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Geneformer compatibility runtime: %s\n' "${FAILED}" >&2
  fi
  if [[ -d "${WORK}" ]]; then
    rm -rf -- "${WORK}"
  fi
  exit "${exit_code}"
}
trap finish EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Geneformer compatibility runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
if [[ ! -x "${SOURCE_ENV}/bin/python" ]]; then
  printf 'Required existing Geneformer Python runtime is absent.\n' >&2
  exit 1
fi
module purge
module load python/3.11.5-GCCcore-13.2.0
CONTROL_PYTHON=$(command -v python)
install -d -m 0750 "${STAGE}/payload" "${STAGE}/validation" "${WORK}/runtime"
cd "${ROOT}"

printf '%s  %s\n' \
  adcce2fd0c9665a5790233e95a170b19b4047c0bbf7218dab21e22bcaa277dbd "${AUDIT}/runtime_audit.json" \
  12fb2742c61e5a02793eb3ea53e84ed4005b1527da894c6b42199753bb3c6d98 "${AUDIT}/pip_freeze.txt" \
  c70e6b43d05a2fc73c8492d97cf0feb716aa4fd6ad42a9c09cf05e3770139654 "${AUDIT}/conda_explicit.txt" \
  b269cef53fa92cc7e7ab9c21a2b85de69904cdfaaf5881e9676c57dc2ae02ae0 "${BUNDLE}/ARTIFACTS.json" \
  b28a589b9e2886776eaf888ea7e82e13c3ef3b25d365a1e4fc572f84a5b4c276 "${INVENTORY}" \
  dece0446a19ce0746e860db615e82a53219c03ebd5cd13ff16580954a6faa38c "${TOKENS}" \
  d7ba1227048487e485bdb80f6b2ee1018f9b9d80b20e9e838beac3f1049930ca "${TOKEN_SUMMARY}" \
  361ac7ec520bd6dcc2a9628af0c45f8b8dd7c12b4cbb895998e4b66c851150aa "${PRIOR_EMBEDDINGS}" \
  dacd7add5a3e30c02f935168905bed26e2a7c326d99efb48d6c89903d5d995c0 "${PRIOR_ROW_ORDER}" \
  08e285ba7cedaf3065d71e8c2537888721c9bb823982ad26ee034512739f0ff8 "${ADAPTER}" \
  a5f2364570129f68e64ae3ab67480ed28372ab52298b29d84b491b91c465ab70 "${EXTRACTOR}" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"

PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - "${BUNDLE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
PY

"${SOURCE_ENV}/bin/python" -m unittest -v \
  tests.unit.test_preflight_geneformer_v2_316m_common \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

tar --sort=name --mtime='UTC 2026-08-24' --owner=0 --group=0 --numeric-owner \
  --pax-option=delete=atime,delete=ctime -C "${SOURCE_ENV}" -cf - . \
  | gzip -1n > "${STAGE}/payload/geneformer-decima-compat-env.tar.gz"
sha256sum "${STAGE}/payload/geneformer-decima-compat-env.tar.gz" \
  > "${STAGE}/validation/runtime-archive.sha256"
tar -xzf "${STAGE}/payload/geneformer-decima-compat-env.tar.gz" -C "${WORK}/runtime"

"${WORK}/runtime/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
"${WORK}/runtime/bin/python" --version \
  > "${STAGE}/validation/effective-python-version.txt" 2>&1
cp "${AUDIT}/runtime_audit.json" "${STAGE}/validation/source-runtime-audit.json"
cp "${AUDIT}/conda_explicit.txt" "${STAGE}/validation/source-conda-explicit.txt"
cp "${AUDIT}/pip_freeze.txt" "${STAGE}/validation/source-pip-freeze.txt"

"${WORK}/runtime/bin/python" \
  "${ROOT}/scripts/preflight_geneformer_v2_316m_common.py" \
  --mode cpu \
  --bundle "${BUNDLE}" \
  --inventory "${INVENTORY}" \
  --tokens "${TOKENS}" \
  --token-summary "${TOKEN_SUMMARY}" \
  --prior-embeddings "${PRIOR_EMBEDDINGS}" \
  --prior-row-order "${PRIOR_ROW_ORDER}" \
  --adapter "${ADAPTER}" \
  --extraction-script "${EXTRACTOR}" \
  --output "${STAGE}/cpu_gate_receipt.json"

sha256sum \
  "${ROOT}/scripts/preflight_geneformer_v2_316m_common.py" \
  "${ROOT}/tests/unit/test_preflight_geneformer_v2_316m_common.py" \
  "${ROOT}/slurm/build_geneformer_v2_316m_compat_runtime_cpu.sbatch" \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "cpu_gate_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_cpu_runtime_and_strict_restore"
    or receipt.get("raw_unrefined_embeddings_only") is not True
    or receipt.get("labels_loaded") is not False
    or receipt.get("sealed_outcomes_loaded") is not False
    or receipt.get("result", {}).get("strict_to_prior_encoder_parameter_identity") is not True
):
    raise SystemExit("Geneformer CPU compatibility receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "geneformer_v2_316m_compatibility_runtime",
        "model_id": "geneformer_v2_316m",
        "runtime_archive_sha256": sha256_file(
            stage / "payload/geneformer-decima-compat-env.tar.gz"
        ),
        "runtime_claim": "compatibility_snapshot_not_upstream_canonical_environment",
        "strict_checkpoint_restore": True,
        "strict_to_prior_encoder_parameter_identity": True,
        "raw_unrefined_embeddings_only": True,
        "prior_smoke_numeric_reproduction_completed": False,
        "project_labels_read": False,
        "sealed_outcomes_read": False,
        "status": "passed_cpu_gate",
    },
)
PY

trap - EXIT
rm -rf -- "${WORK}"
mv "${STAGE}" "${TARGET}"
printf 'Geneformer V2-316M compatibility runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
