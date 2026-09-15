#!/usr/bin/env bash
#SBATCH --job-name=model-data-077
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00
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
TARGET=${ROOT}/executions/model-data-077-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite RSEM-reference artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load RSEM/1.3.3-foss-2022a STAR/2.7.10b-GCC-11.3.0
GSE268273_TOOL_PATH=${PATH}
GSE268273_TOOL_LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
rsem-calculate-expression --version > "${STAGE}/rsem_version.txt" 2>&1
STAR --version > "${STAGE}/star_version.txt" 2>&1
grep -Fx 'Current version: RSEM v1.3.1' "${STAGE}/rsem_version.txt"
grep -Fx '2.7.10b' "${STAGE}/star_version.txt"
module list > "${STAGE}/tool_module_list.txt" 2>&1
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
export GSE268273_TOOL_PATH
export GSE268273_TOOL_LD_LIBRARY_PATH
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
"${PYTHON}" - "${RUNTIME}" > "${STAGE}/pinned_python_runtime.json" <<'PY'
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
bash -n slurm/build_gse268273_rsem_reference_cpu.sbatch
"${PYTHON}" - "${PLAN_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
PY
printf '%s  %s\n' \
  '29566d52311a4204adbfd7ec4501bf07b452ca770ff81f47a3da47fc89b2d717' \
  "${PLAN_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' \
  '2b32cd901ca34440ae5256a0536df34a14bb8b6651ed955ec2b39dfc65df37b7' \
  "${REVISION}" | sha256sum --check --strict
curl --fail --location --retry 8 --retry-all-errors --connect-timeout 30 \
  --output "${STAGE}/gencode.v36.primary_assembly.annotation.gtf.gz" "${SOURCE_GTF_URL}"
"${PYTHON}" scripts/build_gse268273_rsem_reference.py \
  --genome "${GENOME}" \
  --source-gtf "${STAGE}/gencode.v36.primary_assembly.annotation.gtf.gz" \
  --target-gtf "${TARGET_GTF}" \
  --activation-root "${ACTIVATION}" \
  --model-input-root "${MODEL_INPUT}" \
  --work "${STAGE}/work" \
  --output "${STAGE}/artifact" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
module list > "${STAGE}/python_module_list.txt" 2>&1
sha256sum \
  scripts/plan_gse268273_raw_campaign.py \
  scripts/build_gse268273_rsem_reference.py \
  slurm/build_gse268273_rsem_reference_cpu.sbatch \
  "${REVISION}" \
  "${PLAN_ROOT}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  "${MODEL_INPUT}/ARTIFACTS.json" \
  "${STAGE}/artifact/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {"artifact_class": "gse268273_rsem_reference_execution", "status": "passed"})
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/artifact/ARTIFACTS.json"
