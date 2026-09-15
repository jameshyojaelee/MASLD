#!/usr/bin/env bash
#SBATCH --job-name=source-audit
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=executions/source-audit-%j.out
#SBATCH --error=executions/source-audit-%j.err

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
GEO=${ROOT}/executions/geo-multimodal-source-audit-21064310
JOIN=${ROOT}/executions/gse105127-authoritative-join-rights-21065552
RRBS=${ROOT}/executions/gse105127-rrbs-source-qc-r2-21066008
FINAL=${ROOT}/executions/gse105127-assay-native-activation-${SLURM_JOB_ID}
STAGE=${FINAL}.staging
FAILED=${FINAL}.failed
ENA_URL='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP120376&result=read_run&fields=study_accession%2Csample_accession%2Cexperiment_accession%2Crun_accession%2Cscientific_name%2Clibrary_strategy%2Clibrary_source%2Clibrary_selection%2Clibrary_layout%2Cinstrument_platform%2Cinstrument_model%2Cread_count%2Cbase_count%2Cfastq_ftp%2Cfastq_bytes%2Cfastq_md5&format=tsv&download=true'

cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${ROOT}/src:${ROOT}"
if [[ -e "${FINAL}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to replace an existing activation attempt\n' >&2
  exit 1
fi
mkdir -p "${STAGE}/source"
preserve_failed_attempt() {
  status=$?
  if [[ ${status} -ne 0 && ! -e "${FAILED}" ]]; then
    if [[ -d "${STAGE}" ]]; then
      mv "${STAGE}" "${FAILED}"
    elif [[ -d "${FINAL}" ]]; then
      mv "${FINAL}" "${FAILED}"
    fi
  fi
  exit "${status}"
}
trap preserve_failed_attempt EXIT
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

curl --fail --location --retry 3 --output "${STAGE}/source/ena_runinfo.tsv" "${ENA_URL}"
bash -n slurm/build_gse105127_assay_native_activation_io.sbatch
"${PYTHON}" -m unittest \
  tests.unit.test_build_gse105127_assay_native_activation
"${PYTHON}" scripts/build_gse105127_assay_native_activation.py \
  --geo-root "${GEO}" \
  --join-root "${JOIN}" \
  --rrbs-root "${RRBS}" \
  --ena-runinfo "${STAGE}/source/ena_runinfo.tsv" \
  --output "${STAGE}/fixture"
(
  cd "${STAGE}"
  sha256sum source/ena_runinfo.tsv > artifact-relative.sha256
)
sha256sum \
  scripts/build_gse105127_assay_native_activation.py \
  config/evaluation/gse105127_adjacent_section_rna_rrbs_task.toml \
  tests/unit/test_build_gse105127_assay_native_activation.py \
  slurm/build_gse105127_assay_native_activation_io.sbatch \
  executions/geo-multimodal-source-audit-21064310/ARTIFACTS.json \
  executions/gse105127-authoritative-join-rights-21065552/ARTIFACTS.json \
  executions/gse105127-rrbs-source-qc-r2-21066008/ARTIFACTS.json \
  executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json \
  > "${STAGE}/external-stable-paths.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
freeze_tree(
    root,
    {
        "artifact_class": "gse105127_assay_native_activation_execution",
        "status": "design_ready_no_fit",
        "fit_or_score_performed": False,
    },
)
verify_frozen_tree(root)
PY

mv "${STAGE}" "${FINAL}"
(
  cd "${FINAL}"
  sha256sum --check --strict artifact-relative.sha256
)
(
  cd "${ROOT}"
  sha256sum --check --strict "${FINAL}/external-stable-paths.sha256"
)
"${PYTHON}" - "${FINAL}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
PY
trap - EXIT
printf '%s\n' "${FINAL}"
