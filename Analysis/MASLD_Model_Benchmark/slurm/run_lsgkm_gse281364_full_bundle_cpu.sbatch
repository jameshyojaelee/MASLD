#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-train-614
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=5
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

: "${SPLIT_ID:?SPLIT_ID is required}"
: "${SEED_SET:?colon-delimited SEED_SET is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/lsgkm_gse281364_dinucleotide_full_campaign.json"
SCRIPT="${ROOT}/scripts/run_lsgkm_gse281364_fit.py"
AUTH_SCRIPT="${ROOT}/scripts/freeze_lsgkm_gse281364_full_campaign_authorization.py"
SCALE_HELPER="${ROOT}/scripts/run_lsgkm_gse281364_scale_probe.py"
ADAPTER="${ROOT}/scripts/lsgkm_safe_adapter.py"
TEST="${ROOT}/tests/unit/test_run_lsgkm_gse281364_fit.py"
AUTH_TEST="${ROOT}/tests/unit/test_freeze_lsgkm_gse281364_full_campaign_authorization.py"
AUTHORITY="${ROOT}/executions/model-check-613-21099873/ARTIFACTS.json"
AUTHORITY_SHA256=662aae6634c47821a69a0300f9c0aeed13169d04a591b4c0689db99f3048bb3e
TARGET="${ROOT}/executions/${SLURM_JOB_NAME}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite LS-GKM full bundle: %s\n' "${TARGET}" >&2
  exit 1
fi

IFS=: read -r -a SEEDS <<< "${SEED_SET}"
if [[ "${#SEEDS[@]}" -ne "${SLURM_CPUS_PER_TASK}" ]]; then
  printf 'Seed count and allocated CPU count differ\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/fits"

printf '%s  %s\n' 79d7ca90cde98db5152f544ad25075a8dec0d7f421b9ce710c144880db8aa236 "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 951eab2e2f20099709ce88548e2fd92c14308d21e83a425280bd970ea89f4003 "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' ef27670174ef8697c76064e28cea171dda75650ebc106653230b1dae4fa7466c "${AUTH_SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' 57f81394dc4d6b91cbf69809c47c185719cc24521c439680800477109411560f "${SCALE_HELPER}" | sha256sum --check --strict
printf '%s  %s\n' bef248a9f9c084e503b24d4f1dfb91774287a84ffd5cf1afa2b1c1c97e48596f "${ADAPTER}" | sha256sum --check --strict
printf '%s  %s\n' 8728f8153f3c0b5b7f9a054a57eff3751039ac3e74f50be8e74868faaca9e019 "${TEST}" | sha256sum --check --strict
printf '%s  %s\n' d26bf0f0bd5aee6f165f6722d0d56a4d0bcaa775a11e253f39ef577bc423e9e8 "${AUTH_TEST}" | sha256sum --check --strict
printf '%s  %s\n' "${AUTHORITY_SHA256}" "${AUTHORITY}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_freeze_lsgkm_gse281364_full_campaign_authorization \
  tests.unit.test_run_lsgkm_gse281364_fit \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"

pids=()
for seed in "${SEEDS[@]}"; do
  fit="${STAGE}/fits/${SPLIT_ID}.seed${seed}"
  PYTHONPATH="${ROOT}" python "${SCRIPT}" \
    --project-root "${ROOT}" --config "${CONFIG}" \
    --authorization-artifacts "${AUTHORITY}" \
    --authorization-artifacts-sha256 "${AUTHORITY_SHA256}" \
    --split-id "${SPLIT_ID}" --model-seed "${seed}" \
    --output "${fit}" \
    > "${STAGE}/fits/${SPLIT_ID}.seed${seed}.stdout.json" \
    2> "${STAGE}/fits/${SPLIT_ID}.seed${seed}.stderr.txt" &
  pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then failed=1; fi
done
if [[ "${failed}" -ne 0 ]]; then
  printf 'At least one LS-GKM fit failed\n' >&2
  exit 1
fi

python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${AUTH_SCRIPT}" "${SCALE_HELPER}" "${ADAPTER}" \
  "${TEST}" "${AUTH_TEST}" "${CONFIG}" "${AUTHORITY}" \
  "${ROOT}/slurm/run_lsgkm_gse281364_full_bundle_cpu.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${SPLIT_ID}" "${SEED_SET}" "${AUTHORITY_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
split_id = sys.argv[2]
seeds = [int(value) for value in sys.argv[3].split(":")]
authority_sha = sys.argv[4]
fits = []
for seed in seeds:
    fit_root = root / f"fits/{split_id}.seed{seed}"
    receipt = json.loads((fit_root / "fit_receipt.json").read_text(encoding="utf-8"))
    if (
        receipt["status"] != "pass_one_authorized_full_campaign_fit"
        or receipt["split_id"] != split_id
        or receipt["model_seed"] != seed
        or receipt["authorization_artifacts_sha256"] != authority_sha
        or receipt["positive_records"] != 10000
        or receipt["negative_records"] != 10000
        or receipt["held_elements_scored"] <= 0
        or receipt["production_fits_executed"] != 1
        or receipt["benchmark_metrics_calculated"]
        or any(receipt[field] for field in ("outcomes_read", "reporter_counts_read", "control_prediction_values_read", "sealed_assets_read"))
    ):
        raise SystemExit(f"LS-GKM full fit receipt differs: {split_id} seed={seed}")
    fits.append({
        "split_id": split_id,
        "model_seed": seed,
        "model_sha256": receipt["model_receipt"]["sha256"],
        "scores_sha256": receipt["held_element_scores_sha256"],
        "held_elements_scored": receipt["held_elements_scored"],
        "maximum_resident_set_kbytes": receipt["maximum_resident_set_kbytes"],
        "total_native_wall_seconds": receipt["total_native_wall_seconds"],
    })
bundle = {
    "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-fit-bundle-v1",
    "status": "pass_authorized_split_bundle",
    "split_id": split_id,
    "model_seeds": seeds,
    "fit_count": len(fits),
    "fits": fits,
    "authorization_artifacts_sha256": authority_sha,
    "benchmark_metrics_calculated": False,
    "outcomes_read": False,
    "reporter_counts_read": False,
    "control_prediction_values_read": False,
    "sealed_assets_read": False,
}
(root / "bundle_receipt.json").write_text(
    json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exact_dinucleotide_split_fit_bundle",
    "split_id": split_id,
    "model_seeds": seeds,
    "fit_count": len(fits),
    "authorization_artifacts_sha256": authority_sha,
    "benchmark_metrics_calculated": False,
    "outcomes_read": False,
    "reporter_counts_read": False,
    "control_prediction_values_read": False,
    "sealed_assets_read": False,
    "status": bundle["status"],
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM full fit bundle: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
