#!/usr/bin/env bash
#SBATCH --job-name=model-data-880
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
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE_TREE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE135251
COUNTS=${SOURCE_TREE}/GSE135251_counts.tsv
METADATA=${SOURCE_TREE}/GSE135251_metadata.tsv
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
GTF_SHA256=73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9
TARGET=${ROOT}/executions/model-data-880-${SLURM_JOB_ID}-gse135251-source
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE135251 source-freeze artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${GTF_SHA256}" "${GTF}" | sha256sum --check --strict
bash -n slurm/freeze_gse135251_nash_source_inputs_cpu.sbatch
# The source freeze may never name the external series' expression or labels.
if grep -n 'gse83452_gpl16686_gene_matrix\|gse83452_records\|model-data-088\|model-data-072' \
    scripts/freeze_gse135251_nash_source_inputs.py; then
  printf 'The source freeze unexpectedly names a GSE83452 artifact.\n' >&2
  exit 1
fi
"${PYTHON}" -m unittest -v \
  tests.unit.test_freeze_gse135251_nash_source_inputs \
  tests.unit.test_build_gse83452_baseline_selection \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/freeze_gse135251_nash_source_inputs.py \
  --counts "${COUNTS}" \
  --metadata "${METADATA}" \
  --gencode-gtf "${GTF}" \
  --gencode-gtf-sha256 "${GTF_SHA256}" \
  --outer-folds 5 \
  --fold-seed 1701 \
  --output "${STAGE}/source" \
  > "${STAGE}/source_freeze_stdout.json" 2> "${STAGE}/source_freeze_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/freeze_gse135251_nash_source_inputs.py \
  tests/unit/test_freeze_gse135251_nash_source_inputs.py \
  slurm/freeze_gse135251_nash_source_inputs_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${COUNTS}" \
  "${METADATA}" \
  "${GTF}" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
source = root / "source"
receipt = json.loads((source / "source_freeze_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_source_frozen_hashed_copy"
    or receipt.get("participants") != 180
    or receipt.get("metadata_rows") != 180
    or receipt.get("counts_columns_without_metadata") != 36
    or receipt.get("counts_rows_deposited") != 77078
    or receipt.get("external_labels_read") is not False
    or receipt.get("external_expression_values_read") is not False
    or receipt.get("age_sex_logistic_fittable_on_this_source") is not False
    or receipt.get("gate_eligible_arm_id") != "nash_vs_nafl"
    or receipt["arms"]["nash_vs_nafl"]["class_counts"] != {"no_nash": 41, "nash": 131}
    or receipt["arms"]["nash_vs_nafl_plus_control"]["class_counts"]
    != {"no_nash": 49, "nash": 131}
    or receipt["declared_census_reproduced"]["Fibrosis_stage"]
    != {"0": 35, "1": 41, "2": 48, "3": 44, "4": 12}
):
    raise SystemExit("GSE135251 source-freeze verification differs")
for name in ("raw", "molecular", "outcomes", "folds"):
    verify_frozen_tree(source / name)
freeze_tree(
    source,
    {
        "artifact_class": "gse135251_nash_transfer_source",
        "source_series": "GSE135251",
        "participants": 180,
        "source_genes_resolved": receipt["source_genes_resolved"],
        "external_labels_read": False,
        "external_expression_values_read": False,
        "source_roster_admission_state": "proposed_not_approved",
        "status": "passed",
    },
)
verify_frozen_tree(source)
freeze_tree(
    root,
    {
        "artifact_class": "gse135251_nash_transfer_source_campaign",
        "external_labels_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
print(json.dumps({
    "source_genes_resolved": receipt["source_genes_resolved"],
    "counts_row_key_states": receipt["counts_row_key_states"],
    "arms": {k: v["class_counts"] for k, v in receipt["arms"].items()},
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/source/ARTIFACTS.json" \
  "${TARGET}/source/raw/ARTIFACTS.json" \
  "${TARGET}/source/molecular/ARTIFACTS.json" \
  "${TARGET}/source/outcomes/ARTIFACTS.json" \
  "${TARGET}/source/folds/ARTIFACTS.json"
