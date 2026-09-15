#!/usr/bin/env bash
#SBATCH --job-name=model-evaluation
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
EMBEDDINGS="${ROOT}/executions/model-work-301-21072684"
EMBEDDINGS_SHA=641ef6be0436e323a473dfacbd7f66f673cafc498a795bdcddf9b048e44d1b20
OUTCOMES="${ROOT}/executions/gse281364-replicate-outcomes-21066307"
OUTCOMES_SHA=cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db
FIXTURE="${ROOT}/executions/gse281364-caduceus-131k-fixture-21071004"
FIXTURE_SHA=82b47ba2059193cae23e35db5401013659879024052c47558e602c03b4ee2c70
PREFLIGHT="${ROOT}/executions/caduceus-evaluator-preflight-21072511"
PREFLIGHT_SHA=9ac8e136b28143388dd1bb1f17a40e2b6dabab036d52a4c7083c286579f2be5e
TARGET="${ROOT}/executions/gse281364-caduceus-common-lane-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Caduceus evaluation artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for binding in \
  "${EMBEDDINGS_SHA} ${EMBEDDINGS}" \
  "${OUTCOMES_SHA} ${OUTCOMES}" \
  "${FIXTURE_SHA} ${FIXTURE}" \
  "${PREFLIGHT_SHA} ${PREFLIGHT}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" | sha256sum --check --strict
done
for binding in \
  "bb837a6abbf6d01a504d9539345c970d1fec14d37f19b6a7b75431251dd3e034 ${ROOT}/scripts/evaluate_gse281364_caduceus_common_lane.py" \
  "77d9bdfeedde8ad55ca5858fc7325cdb80fe59e1d27d9737abe8eb32b34d1b7c ${ROOT}/tests/unit/test_evaluate_gse281364_caduceus_common_lane.py" \
  "39730c63b7370ea05153482f881da5165287e0a8ff1ad757a1d9399db2f7ab82 ${ROOT}/scripts/evaluate_gse281364_dna_lm_common_lane.py" \
  "066d37b0497008caac7b96006ec2b2ebab132afa06189bb518cab435a1ab9afc ${ROOT}/tests/unit/test_evaluate_gse281364_dna_lm_common_lane.py"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${EMBEDDINGS}" "${OUTCOMES}" "${FIXTURE}" "${PREFLIGHT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
module purge
cd "${ROOT}"
PYTHONPATH="${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse281364_dna_lm_common_lane \
  tests.unit.test_evaluate_gse281364_caduceus_common_lane \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" "${PYTHON}" scripts/evaluate_gse281364_caduceus_common_lane.py \
  --embeddings "${EMBEDDINGS}" \
  --embeddings-sha256 "${EMBEDDINGS_SHA}" \
  --outcomes "${OUTCOMES}" \
  --outcomes-sha256 "${OUTCOMES_SHA}" \
  --fixture "${FIXTURE}" \
  --fixture-sha256 "${FIXTURE_SHA}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json" \
  2> "${STAGE}/validation/evaluation.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  scripts/evaluate_gse281364_caduceus_common_lane.py \
  tests/unit/test_evaluate_gse281364_caduceus_common_lane.py \
  scripts/evaluate_gse281364_dna_lm_common_lane.py \
  tests/unit/test_evaluate_gse281364_dna_lm_common_lane.py \
  slurm/evaluate_gse281364_caduceus_common_lane_cpu.sbatch \
  "${PREFLIGHT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${EMBEDDINGS_SHA}" <<'PY'
from pathlib import Path
import csv
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
embedding_sha = sys.argv[2]
receipt = json.loads((root / "evaluation/receipt.json").read_text())
audit = json.loads((root / "evaluation/input_audit.json").read_text())
with (root / "evaluation/aggregate_metrics.tsv").open(newline="", encoding="utf-8") as handle:
    metrics = list(csv.DictReader(handle, delimiter="\t"))
with (root / "evaluation/promotion.tsv").open(newline="", encoding="utf-8") as handle:
    promotion = list(csv.DictReader(handle, delimiter="\t"))
with (root / "evaluation/head_selection.tsv").open(newline="", encoding="utf-8") as handle:
    selection = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt["status"] != "pass_descriptive_one_seed_development_smoke"
    or receipt["models"] != ["caduceus"]
    or receipt["elements"] != 1033
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["experimental_replicates_per_context"] != 4
    or receipt["experimental_replicates_used_as_independent_donors"]
    or receipt["donor_count"] != 0
    or receipt["sealed_outcomes_read"]
    or receipt["external_evaluation"]
    or receipt["champion_eligible"]
    or receipt["clinical_claim_supported"]
    or receipt["cell_type_eQTL_or_ieQTL_head_fit"]
    or audit["status"] != "pass_independent_projection_consistency_audit"
    or not audit["recursive_member_checksums_verified"]
    or len(audit["checks"]) != 5
    or len(metrics) != 10
    or len(promotion) != 1
    or len(selection) != 30
):
    raise SystemExit("Caduceus common-lane evaluation contract differs")
freeze_tree(root, {
    "artifact_class": "gse281364_caduceus_common_lane_development_evaluation",
    "dataset_id": "gse281364",
    "model_id": "caduceus",
    "embedding_artifacts_sha256": embedding_sha,
    "outcome_artifacts_sha256": "cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db",
    "fixture_artifacts_sha256": "82b47ba2059193cae23e35db5401013659879024052c47558e602c03b4ee2c70",
    "outcome_role": "exposed_development_MPRA_only",
    "experimental_replicates_per_context": 4,
    "replicates_used_as_independent_donors": False,
    "donor_count": 0,
    "seed_count": 1,
    "external_evaluation": False,
    "champion_eligible": False,
    "status": "passed",
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 Caduceus common-lane evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
