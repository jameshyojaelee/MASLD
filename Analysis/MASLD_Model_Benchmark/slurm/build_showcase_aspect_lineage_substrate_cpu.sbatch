#!/usr/bin/env bash
#SBATCH --job-name=model-data-950
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python

ARM_A_FIXTURE=${ROOT}/executions/model-data-064-21079902/fixture
ARM_A_ACTIVATION=${ROOT}/executions/model-data-061-21079623/activation
ARM_B_SOURCE=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
ARM_C_SPEC=${ROOT}/executions/gse202379-pooled-fibrosis-spec-21114717

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ARM_B_SRA=${PROJECT}/Analysis/Deconvolution/metadata/GSE135251_sra_metadata.tsv
ARM_C_SERIES=${PROJECT}/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/geo_cache/GSE202379_series_matrix.txt.gz
ARM_C_RUNINFO=${PROJECT}/data/GSE202379/metadata/SraRunInfo.csv
ARM_C_PAIRING=${PROJECT}/data/GSE202379/metadata/donor_pairing.csv
ARM_C_H5AD=${PROJECT}/Analysis/SingleCell/results_gpu_v2_phase05/atlas/per_dataset/GSE202379_v2.h5ad
BAYESPRISM=${PROJECT}/Analysis/Deconvolution/results
BAYESPRISM_MANIFEST=${ROOT}/config/showcase_substrate/bayesprism_sources.sha256

TARGET=${ROOT}/executions/model-data-950-${SLURM_JOB_ID}-showcase-substrate
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite showcase-substrate artifacts: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
hostname
"${PYTHON}" -V
bash -n slurm/build_showcase_aspect_lineage_substrate_cpu.sbatch

# Every frozen input is checked before a byte of it is read.  A drifted input is
# a failed run, not a quietly different substrate.
cat > "${STAGE}/validation/inputs.sha256" <<'SHA'
0a0604a52907d5bd022ee16128d19657f88fe82080afa309ecc9d8f9638e4002  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/ARTIFACTS.json
e62c882baaaa236af2a6d4a6818abc458ad0810c30014969ccee9d055b8df299  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-data-061-21079623/ARTIFACTS.json
4dab6a759371ccca56227bc3ce1df2c4e9ad27397cab07db4c337669c5003e8d  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-data-880-21130257-gse135251-source/ARTIFACTS.json
e30d8b423b83e325530ebebfb6bce2dca0de8306a04d2a9b66b3151e02b83fc1  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/gse202379-pooled-fibrosis-spec-21114717/ARTIFACTS.json
f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json
2994ba6b7cddad4b5f44626fbc3e12f495d6c617db7b70a72770a75e49f060cc  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/geo_cache/GSE202379_series_matrix.txt.gz
5a6e8973b217a19726b90fcafc33b58280b676684c42e753ace0b640536002b5  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/GSE202379/metadata/SraRunInfo.csv
3ce1f50fe863fab210ed097c46c99ee73abf4acded86f2586c02efcf402918fb  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/GSE202379/metadata/donor_pairing.csv
0e75d3ccf98190e11e5fa89473d85ef6465d4171d7e676aeab1258a49e454291  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/metadata/GSE135251_sra_metadata.tsv
c07eafbcf506e0d38f55ce1cd0f30ee8274b756f74736a628d8ce8edeb26dbff  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2_phase05/atlas/per_dataset/GSE202379_v2.h5ad
SHA
sha256sum --check --strict "${STAGE}/validation/inputs.sha256"
sha256sum --check --strict "${BAYESPRISM_MANIFEST}"
cp "${BAYESPRISM_MANIFEST}" "${STAGE}/validation/bayesprism_sources.sha256"

# The frozen input trees must verify against their own manifests too.
PYTHONPATH="${ROOT}/src" "${PYTHON}" - <<'PY'
from pathlib import Path
from masld_bench.artifacts import verify_frozen_tree

root = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark")
for tree in (
    "executions/model-data-064-21079902",
    "executions/model-data-061-21079623",
    "executions/model-data-880-21130257-gse135251-source",
    "executions/model-data-880-21130257-gse135251-source/source",
    "executions/gse202379-pooled-fibrosis-spec-21114717",
):
    verify_frozen_tree(root / tree)
    print(f"verified {tree}")
PY

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_build_showcase_aspect_lineage_substrate \
  tests.unit.test_freeze_gse135251_nash_source_inputs \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" \
  scripts/build_showcase_aspect_lineage_substrate.py \
  --arm-a-fixture "${ARM_A_FIXTURE}" \
  --arm-a-activation "${ARM_A_ACTIVATION}" \
  --arm-b-source "${ARM_B_SOURCE}" \
  --arm-b-sra-metadata "${ARM_B_SRA}" \
  --arm-c-series-matrix "${ARM_C_SERIES}" \
  --arm-c-run-info "${ARM_C_RUNINFO}" \
  --arm-c-single-nucleus "${ARM_C_H5AD}" \
  --arm-c-donor-pairing "${ARM_C_PAIRING}" \
  --arm-c-frozen-spec "${ARM_C_SPEC}" \
  --bayesprism-results "${BAYESPRISM}" \
  --bayesprism-manifest "${BAYESPRISM_MANIFEST}" \
  --output "${STAGE}/substrate" \
  > "${STAGE}/build_stdout.json" 2> "${STAGE}/build_stderr.txt"

# Nothing outside the benchmark tree may have been written.
sha256sum --check --strict "${STAGE}/validation/inputs.sha256"
sha256sum --check --strict "${BAYESPRISM_MANIFEST}"
printf 'every read input byte-identical before and after\n'

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/build_showcase_aspect_lineage_substrate.py \
  tests/unit/test_build_showcase_aspect_lineage_substrate.py \
  slurm/build_showcase_aspect_lineage_substrate_cpu.sbatch \
  config/showcase_substrate/bayesprism_sources.sha256 \
  scripts/freeze_gse135251_nash_source_inputs.py \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
substrate = root / "substrate"
summary = json.loads((substrate / "substrate_summary.json").read_text(encoding="utf-8"))
arm_a = json.loads((substrate / "arm_a" / "receipt.json").read_text(encoding="utf-8"))
arm_b = json.loads((substrate / "arm_b" / "receipt.json").read_text(encoding="utf-8"))
arm_c = json.loads((substrate / "arm_c" / "receipt.json").read_text(encoding="utf-8"))
lineage = json.loads((substrate / "lineage" / "receipt.json").read_text(encoding="utf-8"))

# The realised units are asserted against what the deposited sources say, in
# distinct units, not in rows.
if (
    arm_a["distinct_participants"] != 99
    or arm_a["aspects"]["lobular_necrosis_is_a_nas_component"] is not False
    or arm_a["aspects"]["nas_sum_reproduced_from_components"] is not True
    or arm_b["distinct_participants"] != 180
    or arm_b["source_rebuilt"] is not False
    or arm_b["counts_column_reconciliation"]["counts_columns_deposited"] != 216
    or arm_b["counts_column_reconciliation"]["columns_without_bulk_metadata"] != 36
    or arm_b["counts_column_reconciliation"]["columns_without_any_metadata"] != 0
    or arm_c["unit_census"]["distinct_donors"] != 47
    or arm_c["unit_census"]["sample_rows_gsm"] != 59
    or arm_c["unit_census"]["saf_graded_donors"] != 40
    or arm_c["is_training_data"] is not False
    or arm_c["stage_source"]["donor_fstage_documented_tsv_used"] is not False
    or lineage["ref_human_levels"] != 16
    or lineage["row_sums_verified_to_one"] is not True
    or lineage["music_used"] is not False
    or summary["model_fitted"] is not False
    or summary["metrics_calculated"] is not False
):
    raise SystemExit("showcase substrate verification differs from the deposited census")

for name in ("arm_a", "arm_b", "arm_c", "lineage", "shared", "sources"):
    freeze_tree(substrate / name, {
        "artifact_class": f"showcase_substrate_{name}",
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    })
    verify_frozen_tree(substrate / name)
freeze_tree(substrate, {
    "artifact_class": "showcase_aspect_lineage_substrate",
    "arm_a_participants": arm_a["distinct_participants"],
    "arm_b_participants": arm_b["distinct_participants"],
    "arm_c_donors": arm_c["unit_census"]["distinct_donors"],
    "arm_c_saf_graded_donors": arm_c["unit_census"]["saf_graded_donors"],
    "shared_gene_axis": summary["shared_gene_axis"]["size"],
    "lineage_cohorts_frozen": lineage["cohorts_frozen"],
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(substrate)
freeze_tree(root, {
    "artifact_class": "showcase_aspect_lineage_substrate_campaign",
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps(summary["realised_units"], indent=2, sort_keys=True))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'substrate: %s\n' "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/substrate/ARTIFACTS.json" \
  "${TARGET}/substrate/arm_a/ARTIFACTS.json" \
  "${TARGET}/substrate/arm_b/ARTIFACTS.json" \
  "${TARGET}/substrate/arm_c/ARTIFACTS.json" \
  "${TARGET}/substrate/lineage/ARTIFACTS.json" \
  "${TARGET}/substrate/shared/ARTIFACTS.json" \
  "${TARGET}/substrate/sources/ARTIFACTS.json"
