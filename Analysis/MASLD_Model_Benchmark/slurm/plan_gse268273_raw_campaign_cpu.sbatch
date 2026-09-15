#!/usr/bin/env bash
#SBATCH --job-name=model-data-075
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
FIXTURE=${ROOT}/executions/model-data-074-21080832/fixture
PREPROCESSING=${FIXTURE}/preprocessing_input
MODEL_INPUT=${FIXTURE}/model_input
ACTIVATION=${ROOT}/executions/gse268273-activation-21066074
ENA_API=https://www.ebi.ac.uk/ena/portal/api/filereport
TARGET=${ROOT}/executions/model-data-075-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite raw campaign-plan artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
"${PYTHON}" - "${RUNTIME}" > "${STAGE}/pinned_python_before_tool_probe.json" <<'PY'
from pathlib import Path
import json
import sys
import numpy

runtime = Path(sys.argv[1]).resolve(strict=True)
numpy_path = Path(numpy.__file__).resolve(strict=True)
if not numpy_path.is_relative_to(runtime):
    raise SystemExit(f"NumPy escaped the pinned runtime: {numpy_path}")
print(json.dumps({"python": sys.version, "numpy": numpy.__version__, "numpy_path": str(numpy_path)}, sort_keys=True))
PY
bash -n slurm/plan_gse268273_raw_campaign_cpu.sbatch
bash -n slurm/download_gse268273_raw_bundles_io.sbatch
bash -n slurm/build_gse268273_rsem_reference_cpu.sbatch
bash -n slurm/quantify_gse268273_rsem_bundles_cpu.sbatch
bash -n slurm/consolidate_gse268273_rsem_counts_cpu.sbatch
bash -n slurm/submit_gse268273_raw_dag.sh
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse268273_raw_campaign \
  > "${STAGE}/unit_tests.txt" 2>&1
curl --fail --location --retry 8 --retry-all-errors --connect-timeout 30 \
  --get "${ENA_API}" \
  --data-urlencode 'accession=PRJNA1116068' \
  --data-urlencode 'result=read_run' \
  --data-urlencode 'fields=run_accession,experiment_accession,library_layout,read_count,base_count,fastq_ftp,fastq_bytes,fastq_md5' \
  --data-urlencode 'format=tsv' \
  --output "${STAGE}/ena_read_metadata.tsv"
"${PYTHON}" scripts/plan_gse268273_raw_campaign.py \
  --preprocessing-root "${PREPROCESSING}" \
  --model-input-root "${MODEL_INPUT}" \
  --activation-root "${ACTIVATION}" \
  --ena-read-metadata "${STAGE}/ena_read_metadata.tsv" \
  --bundles 8 \
  --output "${STAGE}/plan" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
module purge
module load RSEM/1.3.3-foss-2022a STAR/2.7.10b-GCC-11.3.0
rsem-calculate-expression --version > "${STAGE}/rsem_version.txt" 2>&1
STAR --version > "${STAGE}/star_version.txt" 2>&1
grep -Fx 'Current version: RSEM v1.3.1' "${STAGE}/rsem_version.txt"
grep -Fx '2.7.10b' "${STAGE}/star_version.txt"
command -v curl gzip rsem-prepare-reference rsem-calculate-expression STAR \
  > "${STAGE}/runtime_commands.txt"
module list > "${STAGE}/module_list.txt" 2>&1
GSE268273_TOOL_PATH=${PATH}
GSE268273_TOOL_LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
export GSE268273_TOOL_PATH
export GSE268273_TOOL_LD_LIBRARY_PATH
"${PYTHON}" - "${RUNTIME}" > "${STAGE}/pinned_python_after_tool_probe.json" <<'PY'
from pathlib import Path
import json
import sys
import numpy

runtime = Path(sys.argv[1]).resolve(strict=True)
numpy_path = Path(numpy.__file__).resolve(strict=True)
if not numpy_path.is_relative_to(runtime):
    raise SystemExit(f"NumPy escaped the pinned runtime: {numpy_path}")
print(json.dumps({"python": sys.version, "numpy": numpy.__version__, "numpy_path": str(numpy_path)}, sort_keys=True))
PY
cmp "${STAGE}/pinned_python_before_tool_probe.json" \
  "${STAGE}/pinned_python_after_tool_probe.json"
"${PYTHON}" -c 'import numpy; print(numpy.__version__)' \
  > "${STAGE}/numpy_version.txt"
"${PYTHON}" - "${STAGE}/plan" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

plan = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((plan / "receipt.json").read_text(encoding="utf-8"))
with (plan / "bundles.tsv").open(encoding="utf-8", newline="") as handle:
    bundles = list(csv.DictReader(handle, delimiter="\t"))
with (plan / "participants.tsv").open(encoding="utf-8", newline="") as handle:
    participants = list(csv.DictReader(handle, delimiter="\t"))
with (plan / "fastq_files.tsv").open(encoding="utf-8", newline="") as handle:
    files = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt.get("status") != "planned_not_submitted"
    or receipt.get("participants") != 109
    or receipt.get("technical_runs") != 824
    or receipt.get("fastq_files") != 968
    or receipt.get("fastq_bytes") != 512881099727
    or receipt.get("effective_single_end_runs") != 680
    or receipt.get("effective_paired_end_runs") != 144
    or receipt.get("effective_single_end_participants") != 73
    or receipt.get("effective_paired_end_participants") != 36
    or receipt.get("declared_layout_discrepancy_admitted") is not True
    or receipt.get("labels_included") is not False
    or len(bundles) != 8
    or len(participants) != 109
    or len(files) != 968
    or len({(row["row_id"], row["bundle_id"]) for row in participants}) != 109
    or any(int(row["fastq_bytes"]) > 66_000_000_000 for row in bundles)
):
    raise SystemExit("GSE268273 raw campaign plan verification differs")
freeze_tree(
    plan,
    {
        "artifact_class": "gse268273_raw_campaign_plan",
        "participants": 109,
        "bundles": 8,
        "fastq_files": 968,
        "fastq_bytes": 512881099727,
        "status": "planned_not_submitted",
    },
)
verify_frozen_tree(plan)
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/plan_gse268273_raw_campaign.py \
  scripts/download_gse268273_raw_bundle.py \
  scripts/build_gse268273_rsem_reference.py \
  scripts/quantify_gse268273_rsem_bundle.py \
  scripts/consolidate_gse268273_rsem_counts.py \
  tests/unit/test_gse268273_raw_campaign.py \
  slurm/plan_gse268273_raw_campaign_cpu.sbatch \
  slurm/download_gse268273_raw_bundles_io.sbatch \
  slurm/build_gse268273_rsem_reference_cpu.sbatch \
  slurm/quantify_gse268273_rsem_bundles_cpu.sbatch \
  slurm/consolidate_gse268273_rsem_counts_cpu.sbatch \
  slurm/submit_gse268273_raw_dag.sh \
  config/campaigns/gse268273_raw_reprocessing_v1.toml \
  "${PREPROCESSING}/ARTIFACTS.json" \
  "${MODEL_INPUT}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {"artifact_class": "gse268273_raw_campaign_plan_execution", "status": "passed"})
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/plan/ARTIFACTS.json"
