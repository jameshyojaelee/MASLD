#!/usr/bin/env bash
#SBATCH --job-name=model-data-611
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
CONFIG="${ROOT}/config/lsgkm_gse281364_dinucleotide_input_materialization.json"
BUILDER="${ROOT}/scripts/build_lsgkm_gse281364_dinucleotide_inputs.py"
AUDITOR="${ROOT}/scripts/audit_lsgkm_gse281364_dinucleotide_inputs.py"
SHUFFLER="${ROOT}/scripts/freeze_lsgkm_gse281364_dinucleotide_null_readiness.py"
TEST="${ROOT}/tests/unit/test_build_lsgkm_gse281364_dinucleotide_inputs.py"
TARGET="${ROOT}/executions/model-data-611-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM dinucleotide inputs: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 34910250a5f2dd328158d4be2d5940124237bfbe17fccc71a38e2ae69439d571 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 9f9aba90dbef7b3a96604e0b659b0eadd5de1b37401d5e1e9b06e420cd5011bc "${BUILDER}" | sha256sum --check --strict
printf '%s  %s\n' 9947602a0c05cb02be01187b8e912cfda6e8d90b545f55a9eca226bfb91adcc7 "${AUDITOR}" | sha256sum --check --strict
printf '%s  %s\n' 13bee3576e08194ed26a90be477b5c59446b693fdcb6f63f1bbf5490b6c89d8d "${SHUFFLER}" | sha256sum --check --strict
printf '%s  %s\n' 248544a57aafce11827d9036a245affa0a5081da470f84f0daebb3219b8564aa "${TEST}" | sha256sum --check --strict

module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python -m unittest -v \
  tests.unit.test_build_lsgkm_gse281364_dinucleotide_inputs \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python "${BUILDER}" \
  --project-root "${ROOT}" --config "${CONFIG}" \
  --output "${STAGE}/materialization" \
  > "${STAGE}/builder.stdout.json" 2> "${STAGE}/builder.stderr.txt"
singularity exec --cleanenv --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${CORE}/site-packages" --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
singularity --version > "${STAGE}/singularity.version.txt"
module purge

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python "${AUDITOR}" \
  --materialization "${STAGE}/materialization" \
  --output "${STAGE}/independent_emitted_fasta_audit.json" \
  > "${STAGE}/audit.stdout.json" 2> "${STAGE}/audit.stderr.txt"
sha256sum "${BUILDER}" "${AUDITOR}" "${SHUFFLER}" "${TEST}" "${CONFIG}" \
  "${ROOT}/slurm/build_lsgkm_gse281364_dinucleotide_inputs_io_r2.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads((root / "materialization/contract/materialization_contract.json").read_text(encoding="utf-8"))
audit = json.loads((root / "independent_emitted_fasta_audit.json").read_text(encoding="utf-8"))
if (
    contract["status"] != "inputs_materialized_pending_independent_emitted_fasta_audit"
    or audit["status"] != "pass_independent_emitted_fasta_invariant_audit"
    or audit["fit_count"] != 25
    or audit["audited_pair_rows"] != 250000
    or audit["scoring_elements"] != 1033
    or audit["scoring_long_range_blocks"] != 239
    or not audit["dinucleotide_counts_rechecked_from_emitted_FASTA"]
    or contract["production_fits_executed"] != 0
    or any(contract[field] for field in ("outcomes_read", "reporter_counts_read", "prediction_values_read", "sealed_assets_read"))
):
    raise SystemExit("LS-GKM emitted input audit differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exact_dinucleotide_null_inputs",
    "design_id": contract["design_id"],
    "split_count": 5,
    "shared_fit_count": 25,
    "positive_FASTA_count": 5,
    "negative_FASTA_count": 25,
    "audited_pair_rows": 250000,
    "scoring_elements": 1033,
    "scoring_long_range_blocks": 239,
    "emitted_FASTA_invariant_audit_passed": True,
    "production_training_authorized": False,
    "production_fits_executed": 0,
    "production_predictions_generated": 0,
    "outcomes_read": False,
    "reporter_counts_read": False,
    "prediction_values_read": False,
    "sealed_assets_read": False,
    "status": "passed_inputs_frozen_scale_probe_pending",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM dinucleotide inputs: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
