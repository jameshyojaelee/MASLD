#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-access-d0g0-valid
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
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
MANIFEST="${ROOT}/config/evaluation/sequence_accessibility_d0g0_predictions.tsv"
MANIFEST_SHA=59ddf83b1f8206611234a8b28fa143dacbe900efadd2e539291f7eb0012c165d
OBSERVED="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-profiles-21064545"
OBSERVED_SHA=a0cf1a5d62f4d8aaa78af8f5368e2ee7c33de99cea7867eeda90a74c7e335dba
BASELINE="${ROOT}/executions/sequence-training-mean-hepatocyte-donor0-genomic0-21064587"
BASELINE_SHA=76980c7ba990aa5948feedcfd02f13ef3d3f550f702512b7a7fd1ccceca4fc82
TARGET="${ROOT}/executions/sequence-accessibility-donor0-genomic0-valid-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed sequence accessibility evaluation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite sequence accessibility evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${MANIFEST_SHA}" "${MANIFEST}" | sha256sum --check --strict
for record in \
  "${OBSERVED_SHA} ${OBSERVED}/ARTIFACTS.json" \
  "${BASELINE_SHA} ${BASELINE}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${OBSERVED}" "${BASELINE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("evaluator-only inputs frozen-tree verification passed")
PY
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
cd "${ROOT}"
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_evaluate_sequence_accessibility \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="src:." python scripts/evaluate_sequence_accessibility.py \
  --prediction-manifest "${MANIFEST}" \
  --observed-profiles "${OBSERVED}" \
  --observed-artifacts-sha256 "${OBSERVED_SHA}" \
  --training-baseline "${BASELINE}" \
  --baseline-artifacts-sha256 "${BASELINE_SHA}" \
  --evaluation-role valid \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json"
sha256sum \
  scripts/evaluate_sequence_accessibility.py \
  tests/unit/test_evaluate_sequence_accessibility.py \
  slurm/evaluate_sequence_accessibility_valid_donor0_genomic0.sbatch \
  config/evaluation/sequence_accessibility_d0g0_predictions.tsv \
  > "${STAGE}/source.sha256"

micromamba deactivate
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "evaluation/evaluation.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["evaluation_role"] != "valid"
    or receipt["test_outcomes_read"]
    or receipt["champion_claim_allowed"]
    or receipt["external_evaluation"]
    or receipt["promotion_gate_evaluated"]
    or len(receipt["prediction_artifacts"]) != 4
):
    raise SystemExit("sequence accessibility validation contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "sequence_accessibility_valid_development_evaluation",
        "dataset_id": "gse296875",
        "split_id": "donor0_genomic0",
        "lineage_id": "hepatocyte",
        "evaluation_role": "valid",
        "models": sorted(receipt["prediction_artifacts"]),
        "biological_unit": "donor",
        "test_outcomes_read": False,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Sequence accessibility valid evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
