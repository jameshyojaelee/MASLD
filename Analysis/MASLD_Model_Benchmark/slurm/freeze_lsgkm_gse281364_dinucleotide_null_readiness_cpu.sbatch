#!/usr/bin/env bash
#SBATCH --job-name=model-check-609
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/lsgkm_gse281364_dinucleotide_null_readiness.json"
SCRIPT="${ROOT}/scripts/freeze_lsgkm_gse281364_dinucleotide_null_readiness.py"
TEST="${ROOT}/tests/unit/test_freeze_lsgkm_gse281364_dinucleotide_null_readiness.py"
TARGET="${ROOT}/executions/model-check-609-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite readiness artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' d3ca0e3e66768d3613c52913a0e464f1f16070d1842e8e69f30e2b2ffd55ec93 \
  "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 13bee3576e08194ed26a90be477b5c59446b693fdcb6f63f1bbf5490b6c89d8d \
  "${SCRIPT}" | sha256sum --check --strict
printf '%s  %s\n' b400f3a950ab841065348e93291972218f444f62f1c873f77f0b49d2e849ea1b \
  "${TEST}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
python -c \
  'from pathlib import Path; import sys; compile(Path(sys.argv[1]).read_text(encoding="utf-8"), sys.argv[1], "exec")' \
  "${SCRIPT}"
cd "${ROOT}"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_freeze_lsgkm_gse281364_dinucleotide_null_readiness \
  > "${STAGE}/tests.stdout.txt" \
  2> "${STAGE}/tests.stderr.txt"
python "${SCRIPT}" \
  --project-root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/readiness_receipt.json" \
  > "${STAGE}/audit.stdout.json" \
  2> "${STAGE}/audit.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${CONFIG}" "${TEST}" > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "readiness_receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "ready_for_outcome_blind_input_materialization_and_scale_probe_only"
    or receipt["scoring_universe"]["element_count"] != 1033
    or receipt["scoring_universe"]["long_range_block_count"] != 239
    or receipt["negative_design"]["positive_count_per_fit"] != 10000
    or receipt["synthetic_shuffle_fixture"]["five_seed_canonical_shuffle_count"] != 5
    or receipt["campaign_relation"]["failed_campaign_evidence_reused"]
    or receipt["production_fits_executed"] != 0
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["prediction_values_read"]
    or receipt["sealed_assets_read"]
):
    raise SystemExit("LS-GKM exact-dinucleotide-null readiness differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_exact_dinucleotide_null_readiness",
    "design_id": receipt["design_id"],
    "element_count": 1033,
    "long_range_block_count": 239,
    "fixed_seed_count": 5,
    "training_source_count": 5,
    "positive_count_per_fit": 10000,
    "negative_count_per_fit": 10000,
    "failed_campaign_evidence_reused": False,
    "production_fits_executed": 0,
    "outcomes_read": False,
    "reporter_counts_read": False,
    "prediction_values_read": False,
    "sealed_assets_read": False,
    "status": receipt["status"],
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'LS-GKM exact-dinucleotide-null readiness: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
