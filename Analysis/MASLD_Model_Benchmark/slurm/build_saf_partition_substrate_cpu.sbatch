#!/usr/bin/env bash
#SBATCH --job-name=model-data-951
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
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

SUBSTRATE=${ROOT}/executions/model-data-950-21181338-showcase-substrate/substrate
SUBSTRATE_SHA=8fed068f5c53d8b0ebbba8c98a734ecd9c99253104a21631e9f58b1c3193bb4b
SERIES=${PROJECT}/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/geo_cache/GSE202379_series_matrix.txt.gz
GRADED=${ROOT}/executions/gse202379-pooled-fibrosis-spec-21114717/spec/graded_donors.tsv
PAIRING=${PROJECT}/data/GSE202379/metadata/donor_pairing.csv

TARGET=${ROOT}/executions/model-data-951-${SLURM_JOB_ID}-saf-partition
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
hostname; "${PYTHON}" -V
bash -n slurm/build_saf_partition_substrate_cpu.sbatch

cat > "${STAGE}/validation/inputs.sha256" <<'SHA'
8fed068f5c53d8b0ebbba8c98a734ecd9c99253104a21631e9f58b1c3193bb4b  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-data-950-21181338-showcase-substrate/substrate/ARTIFACTS.json
e30d8b423b83e325530ebebfb6bce2dca0de8306a04d2a9b66b3151e02b83fc1  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/gse202379-pooled-fibrosis-spec-21114717/ARTIFACTS.json
f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json
2994ba6b7cddad4b5f44626fbc3e12f495d6c617db7b70a72770a75e49f060cc  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/geo_cache/GSE202379_series_matrix.txt.gz
3ce1f50fe863fab210ed097c46c99ee73abf4acded86f2586c02efcf402918fb  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/GSE202379/metadata/donor_pairing.csv
SHA
sha256sum --check --strict "${STAGE}/validation/inputs.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - <<'PY'
from pathlib import Path
from masld_bench.artifacts import verify_frozen_tree
root = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark")
for tree in (
    "executions/model-data-950-21181338-showcase-substrate",
    "executions/model-data-950-21181338-showcase-substrate/substrate",
    "executions/model-data-950-21181338-showcase-substrate/substrate/arm_a",
    "executions/model-data-950-21181338-showcase-substrate/substrate/arm_c",
    "executions/gse202379-pooled-fibrosis-spec-21114717",
):
    verify_frozen_tree(root / tree)
    print(f"verified {tree}")
PY

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_build_saf_partition_substrate \
  tests.unit.test_build_showcase_aspect_lineage_substrate \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" \
  scripts/build_saf_partition_substrate.py \
  --substrate "${SUBSTRATE}" \
  --substrate-artifacts-sha256 "${SUBSTRATE_SHA}" \
  --series-matrix "${SERIES}" \
  --graded-donors "${GRADED}" \
  --donor-pairing "${PAIRING}" \
  --output "${STAGE}/saf_partition" \
  > "${STAGE}/build_stdout.json" 2> "${STAGE}/build_stderr.txt"

sha256sum --check --strict "${STAGE}/validation/inputs.sha256"
printf 'every read input byte-identical before and after\n'

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/build_saf_partition_substrate.py \
  scripts/build_showcase_aspect_lineage_substrate.py \
  tests/unit/test_build_saf_partition_substrate.py \
  slurm/build_saf_partition_substrate_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
saf = root / "saf_partition"
r = json.loads((saf / "receipt.json").read_text(encoding="utf-8"))
c = r["arm_c"]["graded_donors_tsv_characterisation"]
if (
    r["model_fitted"] is not False
    or r["metrics_calculated"] is not False
    or r["arm_a"]["distinct_participants"] != 99
    or r["arm_a"]["saf_triple_state"] != "derived"
    or r["arm_a"]["derivation"]["check_passed"] is not True
    or r["arm_a"]["levels_realised"]["activity_A"] != 5
    or r["arm_c"]["saf_graded_donors"] != 40
    or r["arm_c"]["usable_for_expression_join"] != 39
    or r["arm_c"]["saf_triple_state"] != "observed"
    or c["saf_axes_discarded"] != ["S", "A"]
    or c["fibrosis_disagreements"] != {}
    or c["fibrosis_donors_compared"] != 40
    or r["p70"]["donor_id"] != "P70"
    or r["p70"]["consequence"]["realised_n_for_any_expression_join"] != 39
    or any(o["independent"] for o in r["p70"]["omissions"]
           if o["artifact"] == "donor_fstage_documented.tsv")
):
    raise SystemExit("SAF partition verification differs")
freeze_tree(saf, {
    "artifact_class": "saf_partition_substrate",
    "arm_a_participants": 99,
    "arm_c_usable_donors": 39,
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(saf)
freeze_tree(root, {
    "artifact_class": "saf_partition_substrate_campaign",
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps(r["power_statement"], indent=2, sort_keys=True))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'saf partition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/saf_partition/ARTIFACTS.json"
