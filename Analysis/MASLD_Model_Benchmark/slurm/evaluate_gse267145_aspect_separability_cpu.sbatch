#!/usr/bin/env bash
#SBATCH --job-name=model-check-1044
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0 of the MASLD showcase model, job 2 of 2: the aspect-separability check.
#
# The prespecification frozen by model-check-1043 (job 21181265) is verified by
# digest BEFORE any correlation is computed. The pinned digests below were read
# from that completed tree, so the ordering is provable from the output files and
# does not rest on a claim in a report.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/evaluate_gse267145_aspect_separability.py

PRESPEC_ROOT=${ROOT}/executions/model-check-1043-21181265-aspect-separability-prespec
PRESPEC=${PRESPEC_ROOT}/prespecification
ENDPOINTS=${ROOT}/executions/model-data-061-21079623/activation/participant_endpoints.tsv
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv

TARGET=${ROOT}/executions/model-check-1044-${SLURM_JOB_ID}-aspect-separability
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the separability result.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/evaluate_gse267145_aspect_separability_cpu.sbatch

# ---------------------------------------------------------------------------
# The prespecification came first. Prove it before computing anything.
# ---------------------------------------------------------------------------
cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
b6883741067c3245cc8968b6e878b9b2d4eb6377aeede1c193af5bfaa00d9035  executions/model-check-1043-21181265-aspect-separability-prespec/ARTIFACTS.json
5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76  executions/model-check-1043-21181265-aspect-separability-prespec/prespecification/ARTIFACTS.json
4beecaf18ca0a2c5f15f0037e902f2ebe6863e311aa94adc4f475f3ae2dc996a  executions/model-check-1043-21181265-aspect-separability-prespec/prespecification/aspect_separability_prespecification.json
9a13c0fd65db0483566b7495def11ade3339861de091ea29fcbdc4c34f47b5b9  executions/model-data-061-21079623/activation/participant_endpoints.tsv
812d5fce6e7cc88f6939ae493db7b7e7470057c21cf3b219abc39f587ff64ddc  executions/model-data-064-21079902/fixture/molecular/rna_values.npy
baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd  executions/model-data-064-21079902/fixture/molecular/rna_feature_axis.tsv
7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc  executions/model-data-064-21079902/fixture/molecular/participant_axis.tsv
9ddfaf9c58b58a5b239eebbbc19b214210543fb4f60f3f5a7c563a70fbebb0d1  executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

# The prespecification tree must also be internally consistent, not merely
# match a digest this wrapper carries.
"${PYTHON}" - "${PRESPEC}" "${PRESPEC_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for argument in sys.argv[1:]:
    verify_frozen_tree(Path(argument).resolve(strict=True))
print("prespecification tree verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_evaluate_gse267145_aspect_separability \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

find "${MOLECULAR}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_before.txt"

# The script owns its --output directory and refuses to overwrite one.
"${PYTHON}" "${SCRIPT}" \
  --prespecification "${PRESPEC}" \
  --endpoints "${ENDPOINTS}" \
  --molecular "${MOLECULAR}" \
  --folds "${FOLDS}" \
  --output "${STAGE}/separability" \
  --permutations 20000 \
  --permutation-block 250 \
  --label-draws 20000 \
  --splits 200 \
  --c3-floor-draws 2000 \
  --bootstrap 10000 \
  --c3-bootstrap 2000 \
  --seed 20260827 \
  > "${STAGE}/separability_stdout.txt" 2> "${STAGE}/separability_stderr.txt"

find "${MOLECULAR}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_after.txt"
cmp -s "${STAGE}/frozen_before.txt" "${STAGE}/frozen_after.txt" || { printf 'A frozen artifact changed.\n' >&2; exit 1; }
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" \
  slurm/evaluate_gse267145_aspect_separability_cpu.sbatch \
  tests/unit/test_evaluate_gse267145_aspect_separability.py \
  "${PRESPEC}/aspect_separability_prespecification.json" \
  "${ENDPOINTS}" "${FOLDS}" \
  "${MOLECULAR}/rna_values.npy" "${MOLECULAR}/rna_feature_axis.tsv" \
  "${MOLECULAR}/participant_axis.tsv" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
result = root / "separability"
payload = json.loads((result / "aspect_separability.json").read_text(encoding="utf-8"))

if payload["prespec_id"] != "gse267145_aspect_separability_v1":
    raise SystemExit("the result is not against the frozen prespecification")
if payload["decision"] not in {"GO", "STOP"}:
    raise SystemExit("the decision is neither GO nor STOP")
if payload["a_stop_is_a_result_not_a_failure"] is not True:
    raise SystemExit("a stop is not recorded as a result")
if payload["cohort"]["distinct_participants"] != 99:
    raise SystemExit("the participant count is not 99")
if payload["nas_component_identity"]["component_sum_equals_three_aspects"] != 99:
    raise SystemExit("the NAS component identity does not hold")
if payload["nas_component_identity"]["lobular_necrosis_is_a_nas_component"] is not False:
    raise SystemExit("necrosis reached the three aspects")
if payload["controls"]["fibrosis_reported_as_a_fourth_axis_never_folded_in"] is not True:
    raise SystemExit("fibrosis was folded into the three aspects")
criteria = payload["criteria"]
decided = (
    criteria["c1_labels_are_not_redundant"]["passed"]
    and criteria["c2_effective_dimensionality"]["passed"]
    and criteria["c3_aspects_rank_genes_differently"]["passed"]
)
if decided != (payload["decision"] == "GO"):
    raise SystemExit("the decision does not follow from the three criteria")
for entry in payload["reuse_validation"]["association_agreement_with_masld_bench"]:
    if entry["max_absolute_difference_vs_masld_bench_spearman"] > 1e-10:
        raise SystemExit("the vectorized association disagrees with the package")

verify_frozen_tree(result)
freeze_tree(root, {
    "artifact_class": "gse267145_aspect_separability_campaign",
    "stage": "stage_0_gate",
    "job": "2_of_2",
    "decision": payload["decision"],
    "prespecification_verified_before_computation": True,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "decision": payload["decision"],
    "finding": payload["finding"],
    "c1_passed": criteria["c1_labels_are_not_redundant"]["passed"],
    "c2_passed": criteria["c2_effective_dimensionality"]["passed"],
    "c2_participation_ratio": criteria["c2_effective_dimensionality"][
        "participation_ratio"],
    "c3_passed": criteria["c3_aspects_rank_genes_differently"]["passed"],
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/separability/ARTIFACTS.json"
