#!/usr/bin/env bash
#SBATCH --job-name=masld-scglue-retrieval-prepare
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE_ROOT="${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset"
SOURCE_ARTIFACTS_SHA256=db9e1cb847ec83a61ec79f7a206485db11e323f6b1b54fa4e10bbe4f9728aa3a
SOURCE="${SOURCE_ROOT}/gse296875_rna_atac_smoke_1000.h5"
SOURCE_SHA256=a1c05a432cec6894acae20a129ecceb153f191d9911843c86a90e13cd7f8ce14
PREPARED="${ROOT}/executions/multivi-smoke-prepared-21065093"
PREPARED_ARTIFACTS_SHA256=6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
GTF_SHA256=73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9
TARGET="${ROOT}/executions/scglue-retrieval-prepared-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scGLUE retrieval preparation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scGLUE retrieval preparation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${SOURCE_ARTIFACTS_SHA256} ${SOURCE_ROOT}/ARTIFACTS.json" \
  "${SOURCE_SHA256} ${SOURCE}" \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${RUNTIME_ARTIFACTS_SHA256} ${RUNTIME}/ARTIFACTS.json" \
  "${GTF_SHA256} ${GTF}"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
gzip -t "${GTF}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE_ROOT}" "${PREPARED}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${SITE_PACKAGES}" python "${ROOT}/scripts/scglue_prepare_retrieval_smoke.py" \
  --source "${SOURCE}" \
  --prepared "${PREPARED}/prepared" \
  --gtf "${GTF}" \
  --output "${STAGE}/prepared" \
  > "${STAGE}/validation/prepare.stdout.json"
sha256sum \
  "${ROOT}/scripts/scglue_prepare_retrieval_smoke.py" \
  "${ROOT}/slurm/scglue_prepare_retrieval_smoke.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "prepared/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["pairing_topology"] != "same_nucleus"
    or receipt["outer_unit"] != "donor"
    or set(receipt["folds"]) != set(map(str, range(5)))
    or receipt["hidden_pair_map_available_to_model"]
    or not receipt["query_modality_ids_disjoint"]
    or not receipt["query_atac_row_order_permuted"]
    or receipt["outcomes_read"]
    or receipt["retrieval_metrics_calculated"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("scGLUE retrieval preparation contract differs")
freeze_tree(root, {
    "artifact_class": "scglue_retrieval_prepared_folds",
    "dataset_id": "gse296875",
    "pairing_topology": "same_nucleus",
    "outer_unit": "donor",
    "source_artifacts_sha256": "db9e1cb847ec83a61ec79f7a206485db11e323f6b1b54fa4e10bbe4f9728aa3a",
    "prepared_artifacts_sha256": "6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa",
    "reference_gtf_sha256": "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
    "hidden_pair_map_available_to_model": False,
    "retrieval_metrics_calculated": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scGLUE retrieval preparation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
