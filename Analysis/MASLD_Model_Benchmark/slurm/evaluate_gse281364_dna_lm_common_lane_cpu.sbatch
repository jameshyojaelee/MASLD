#!/usr/bin/env bash
#SBATCH --job-name=masld-dna-lm-common-heads
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
EMBEDDINGS="${ROOT}/executions/model-work-001-21070132"
EMBEDDINGS_SHA=e0efbb81731c1bdb2a1f5546ae3e0678eda28cd5960cbea48afce92f286e62da
OUTCOMES="${ROOT}/executions/gse281364-replicate-outcomes-21066307"
OUTCOMES_SHA=cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db
FIXTURE="${ROOT}/executions/gse281364-dna-lm-common-fixture-21069076"
FIXTURE_SHA=74a3a5c36bd1416daf3fe4b63fe317cfa7e2d008b96e2503d4ca43770b4357de
TARGET="${ROOT}/executions/gse281364-dna-lm-common-lane-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed DNA-LM common-lane stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite DNA-LM evaluation artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

PYTHONPATH="${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse281364_dna_lm_common_lane \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

PYTHONPATH="${ROOT}" "${PYTHON}" \
  "${ROOT}/scripts/evaluate_gse281364_dna_lm_common_lane.py" \
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
  "${ROOT}/scripts/evaluate_gse281364_dna_lm_common_lane.py" \
  "${ROOT}/tests/unit/test_evaluate_gse281364_dna_lm_common_lane.py" \
  "${ROOT}/slurm/evaluate_gse281364_dna_lm_common_lane_cpu.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
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
    or receipt["elements"] != 1033
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["experimental_replicates_per_context"] != 4
    or receipt["experimental_replicates_used_as_independent_donors"]
    or receipt["donor_count"] != 0
    or receipt["sealed_outcomes_read"]
    or receipt["external_evaluation"]
    or receipt["champion_eligible"]
    or receipt["clinical_claim_supported"]
    or audit["status"] != "pass_independent_projection_consistency_audit"
    or not audit["recursive_member_checksums_verified"]
    or len(audit["checks"]) != 15
    or len(metrics) != 18
    or len(promotion) != 3
    or len(selection) != 70
):
    raise SystemExit("DNA-LM common-lane evaluation contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_dna_lm_common_lane_development_evaluation",
        "dataset_id": "gse281364",
        "models": ["dnabert2", "hyenadna", "nucleotide_transformer"],
        "embedding_artifacts_sha256": "e0efbb81731c1bdb2a1f5546ae3e0678eda28cd5960cbea48afce92f286e62da",
        "outcome_artifacts_sha256": "cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db",
        "fixture_artifacts_sha256": "74a3a5c36bd1416daf3fe4b63fe317cfa7e2d008b96e2503d4ca43770b4357de",
        "outcome_role": "exposed_development_MPRA_only",
        "experimental_replicates_per_context": 4,
        "replicates_used_as_independent_donors": False,
        "donor_count": 0,
        "seed_count": 1,
        "restricted_comparator": ["nucleotide_transformer"],
        "external_evaluation": False,
        "champion_eligible": False,
        "status": "passed",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 DNA-LM common-lane evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
