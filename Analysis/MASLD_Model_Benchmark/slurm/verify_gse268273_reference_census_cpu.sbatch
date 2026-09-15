#!/usr/bin/env bash
#SBATCH --job-name=model-data-080
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
: "${PLAN_ROOT:?PLAN_ROOT must name the frozen model-data-075 plan directory}"

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
GENOME=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz
TARGET_GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
ACTIVATION=${ROOT}/executions/gse268273-activation-21066074
MODEL_INPUT=${ROOT}/executions/model-data-074-21080832/fixture/model_input
REVISION=${ROOT}/config/campaigns/gse268273_reference_primary_assembly_revision_v2.toml
SOURCE_GTF_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.primary_assembly.annotation.gtf.gz
TARGET=${ROOT}/executions/model-data-080-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite reference-census verification artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/verify_gse268273_reference_census_cpu.sbatch
bash -n slurm/build_gse268273_rsem_reference_cpu.sbatch
printf '%s  %s\n' \
  '29566d52311a4204adbfd7ec4501bf07b452ca770ff81f47a3da47fc89b2d717' \
  "${PLAN_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' \
  '2b32cd901ca34440ae5256a0536df34a14bb8b6651ed955ec2b39dfc65df37b7' \
  "${REVISION}" | sha256sum --check --strict
"${PYTHON}" - "${PLAN_ROOT}" "${REVISION}" "${RUNTIME}" > "${STAGE}/revision_runtime.json" <<'PY'
from pathlib import Path
import json
import sys
import tomllib
import numpy
from masld_bench.artifacts import verify_frozen_tree

plan = Path(sys.argv[1]).resolve(strict=True)
revision_path = Path(sys.argv[2]).resolve(strict=True)
runtime = Path(sys.argv[3]).resolve(strict=True)
verify_frozen_tree(plan)
with revision_path.open("rb") as handle:
    revision = tomllib.load(handle)
expected = {
    "status": "prepared_not_submitted",
    "failed_reference_job": 21082744,
    "source_gtf_bytes": 44_517_553,
    "source_gtf_md5": "5038acd7158686c69a1f88a5e26a5934",
    "source_crosswalk_features": 14_149,
    "target_stable_genes": 14_078,
}
if any(revision.get(key) != value for key, value in expected.items()):
    raise SystemExit("reference revision contract differs")
numpy_path = Path(numpy.__file__).resolve(strict=True)
if not numpy_path.is_relative_to(runtime):
    raise SystemExit(f"NumPy escaped the pinned runtime: {numpy_path}")
print(json.dumps({
    "numpy": numpy.__version__,
    "numpy_path": str(numpy_path),
    "revision_id": revision["revision_id"],
    "python": sys.version,
}, sort_keys=True))
PY
"${PYTHON}" -m unittest tests.unit.test_gse268273_raw_campaign \
  > "${STAGE}/unit_tests.txt" 2>&1
curl --fail --location --retry 8 --retry-all-errors --connect-timeout 30 \
  --output "${STAGE}/gencode.v36.primary_assembly.annotation.gtf.gz" \
  "${SOURCE_GTF_URL}"
gzip -t "${STAGE}/gencode.v36.primary_assembly.annotation.gtf.gz"
printf '%s  %s\n' '5038acd7158686c69a1f88a5e26a5934' \
  "${STAGE}/gencode.v36.primary_assembly.annotation.gtf.gz" \
  | md5sum --check --strict
[[ $(stat -c '%s' "${STAGE}/gencode.v36.primary_assembly.annotation.gtf.gz") == 44517553 ]]
"${PYTHON}" scripts/build_gse268273_rsem_reference.py \
  --genome "${GENOME}" \
  --source-gtf "${STAGE}/gencode.v36.primary_assembly.annotation.gtf.gz" \
  --target-gtf "${TARGET_GTF}" \
  --activation-root "${ACTIVATION}" \
  --model-input-root "${MODEL_INPUT}" \
  --admission-only > "${STAGE}/admission_receipt.json"
"${PYTHON}" - "${STAGE}/admission_receipt.json" <<'PY'
from pathlib import Path
import json
import sys

receipt = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
expected = {
    "status": "passed",
    "source_annotation_gene_ids": 60_719,
    "source_gtf_contigs": 47,
    "source_crosswalk_features_admitted": 14_149,
    "source_primary_scaffold_features_admitted": 7,
    "target_stable_genes_admitted": 14_078,
    "target_duplicate_sum_groups": 11,
    "source_features_masked_absent_v49": 60,
}
if any(receipt.get(key) != value for key, value in expected.items()):
    raise SystemExit("reference admission receipt differs")
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
module list > "${STAGE}/module_list.txt" 2>&1
sha256sum \
  scripts/plan_gse268273_raw_campaign.py \
  scripts/build_gse268273_rsem_reference.py \
  tests/unit/test_gse268273_raw_campaign.py \
  slurm/verify_gse268273_reference_census_cpu.sbatch \
  slurm/build_gse268273_rsem_reference_cpu.sbatch \
  "${REVISION}" \
  "${PLAN_ROOT}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  "${MODEL_INPUT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "gse268273_reference_census_revision_verification",
    "labels_accessed": False,
    "fit_or_score_performed": False,
    "status": "passed",
})
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
