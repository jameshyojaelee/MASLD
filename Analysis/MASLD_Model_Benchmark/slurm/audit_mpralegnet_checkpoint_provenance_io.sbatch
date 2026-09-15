#!/usr/bin/env bash
#SBATCH --job-name=model-check-226
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
CONFIG="${ROOT}/config/gse281364_mpralegnet_provenance_audit.toml"
TARGET="${ROOT}/executions/model-check-226-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed provenance audit: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite provenance audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/downloads" "${STAGE}/validation"

for binding in \
  "6f77ec23b4990efe416a3f69fcacb9ccce942b2fb062f17855fb9bd32e26ca7f executions/mpralegnet-admission-21066221" \
  "33c3ea9fe134157b9f8f29396cd9d0d565243a9f6e44db3f458cd970c486d32c executions/gse281364-mpralegnet-fixture-21069059" \
  "a59c82c5093a61025f39d5fb12c61820711f58e001ee9bc7858bd505140ec5f9 executions/model-work-102-21069877" \
  "fdc1db620152757f74c07904487a62381fcd182a532669b93b0aa1ca5cdfc2ac executions/model-check-225-21091492"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}/ARTIFACTS.json" \
    | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
for relative in (
    "executions/mpralegnet-admission-21066221",
    "executions/gse281364-mpralegnet-fixture-21069059",
    "executions/model-work-102-21069877",
    "executions/model-check-225-21091492",
):
    verify_frozen_tree(root / relative)
print("MPRALegNet provenance authorities verified")
PY
module purge

module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}/src:${ROOT}"
python -m unittest -v tests.unit.test_audit_mpralegnet_checkpoint_provenance \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

curl_flags=(--fail --location --retry 10 --retry-all-errors --connect-timeout 30 \
  --continue-at - --proto '=https' --tlsv1.2)
curl "${curl_flags[@]}" \
  --output "${STAGE}/downloads/zenodo_record.json" \
  "https://zenodo.org/api/records/10558183"
curl "${curl_flags[@]}" \
  --output "${STAGE}/downloads/final_dump.zip" \
  "https://zenodo.org/api/records/10558183/files/final_dump.zip/content"
curl "${curl_flags[@]}" \
  --output "${STAGE}/downloads/human_legnet-main.zip" \
  "https://zenodo.org/api/records/10558183/files/human_legnet-main.zip/content"

python "${ROOT}/scripts/audit_mpralegnet_checkpoint_provenance.py" \
  --project-root "${ROOT}" \
  --config "${CONFIG}" \
  --zenodo-record "${STAGE}/downloads/zenodo_record.json" \
  --final-dump "${STAGE}/downloads/final_dump.zip" \
  --code-archive "${STAGE}/downloads/human_legnet-main.zip" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
rm -f \
  "${STAGE}/downloads/final_dump.zip" \
  "${STAGE}/downloads/human_legnet-main.zip" \
  "${STAGE}/downloads/zenodo_record.json"
rmdir "${STAGE}/downloads"
sha256sum \
  "${ROOT}/config/gse281364_mpralegnet_provenance_audit.toml" \
  "${ROOT}/scripts/audit_mpralegnet_checkpoint_provenance.py" \
  "${ROOT}/tests/unit/test_audit_mpralegnet_checkpoint_provenance.py" \
  "${ROOT}/slurm/audit_mpralegnet_checkpoint_provenance_io.sbatch" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/provenance/receipt.json").read_text())
passed = receipt["provenance_equivalence_pass"]
contract = root / "audit/contract"
if (
    receipt["status"]
    not in {
        "pass_author_zenodo_tensor_exact_provenance",
        "blocked_checkpoint_tensor_provenance_not_equivalent",
    }
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_assets_read"]
    or receipt["metrics_calculated"]
    or receipt["model_fit"]
    or receipt["predictions_generated"]
    or receipt["schema_seed_repeats_are_independent"]
    or receipt["independent_fitted_predictions"] != 1
    or contract.is_dir() != passed
):
    raise SystemExit("provenance terminal receipt differs")
if passed:
    task = json.loads((contract / "pre_scoring_task_spec.json").read_text())
    evaluator = json.loads((contract / "evaluator_contract.json").read_text())
    if (
        task["claims"]["champion_eligible"]
        or task["claims"]["complementarity_eligible"]
        or task["prediction_universe"]["independent_fitted_predictions"] != 1
        or evaluator["execution"]["authorized_in_this_audit"]
        or evaluator["execution"]["outcomes_read"]
        or evaluator["execution"]["metrics_calculated"]
    ):
        raise SystemExit("pre-scoring contract firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "mpralegnet_checkpoint_provenance_equivalence_audit",
        "dataset_id": "gse281364",
        "model_id": "mpralegnet_hepg2_test1_val2",
        "provenance_equivalence_pass": passed,
        "descriptive_native_output_scoring_eligible": receipt[
            "descriptive_native_output_scoring_eligible"
        ],
        "outcomes_read": False,
        "metrics_calculated": False,
        "model_fit": False,
        "predictions_generated": False,
        "standalone_champion_eligible": False,
        "status": "passed" if passed else "terminal_blocked",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'MPRALegNet provenance audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
