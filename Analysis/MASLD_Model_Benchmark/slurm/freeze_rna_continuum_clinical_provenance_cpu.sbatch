#!/usr/bin/env bash
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
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
COHORT_ROOT="${PROJECT}/RNA-seq/Human/Patient_Cohorts"

# The substrate manifest is byte-identical across all ten figures/candidates/* copies,
# which is what makes this one canonical rather than one copy among ten. Pinning the
# digest also protects this freeze from the fig4f/fig4h adoption writing into that
# tree concurrently.
MANIFEST="${PROJECT}/figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis/stage_extensions/five_cohort_sample_manifest.tsv"
MANIFEST_SHA256=7a446702d0f8995af9f2f90b02f38bb23f3a033bef8312fea51212825caedc0f

TARGET="${ROOT}/executions/rna-continuum-clinical-provenance-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite clinical provenance artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
mkdir -p "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"

# Verify the read-only clinical sources before copying anything. The two archived files
# are unchecksummed at source and are the reason this job exists.
cat > "${STAGE}/validation/bound_fixtures.sha256" <<SUMS
ddf6680858a1c3d0cb25783c51b77bca42b07e19528214b2b02e68e78636f3e4  ${COHORT_ROOT}/archive/old_data_metadata/metadata/GSE130970_SraRunTable.csv
c7488ce83eae2d6abca2e75929d5a3c8550f84af0da8e914562287ae40a44958  ${COHORT_ROOT}/archive/old_data_metadata/metadata/GSE135251_SraRunTable.csv
${MANIFEST_SHA256}  ${MANIFEST}
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

python3 -m unittest -v tests.unit.test_freeze_rna_continuum_clinical_provenance \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

python3 scripts/freeze_rna_continuum_clinical_provenance.py \
  --cohort-root "${COHORT_ROOT}" \
  --dataset-config "${PROJECT}/config/human_datasets.yaml" \
  --harmonize-script "${COHORT_ROOT}/analysis/integration/scripts/00_harmonize_metadata.R" \
  --manifest "${MANIFEST}" \
  --manifest-sha256 "${MANIFEST_SHA256}" \
  --output "${STAGE}" \
  > "${STAGE}/validation/reproduction_stdout.json"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum scripts/freeze_rna_continuum_clinical_provenance.py \
          tests/unit/test_freeze_rna_continuum_clinical_provenance.py \
          slurm/freeze_rna_continuum_clinical_provenance_cpu.sbatch > "${STAGE}/source.sha256"

python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "provenance" / "audit.json").read_text())

if audit["donor_level_claim_made"] is not False:
    raise SystemExit("this artifact may not make a donor-level claim")
if audit["unit_of_observation"] != "sequencing_sample":
    raise SystemExit("the unit of observation must be named as the sequencing sample")
if audit["n_disagreements"] != 0:
    raise SystemExit("the manifest clinical join did not reproduce")
if audit["n_reproduced"] != audit["n_manifest_samples"]:
    raise SystemExit("not every manifest sample was re-derived")

freeze_tree(root, {
    "artifact_class": "rna_continuum_clinical_provenance_freeze",
    "workstream": "H_continuum_rna_descriptive",
    "unit_of_observation": "sequencing_sample",
    "donor_level_claim_made": False,
    "n_donors_status": audit["n_donors_status"],
    "join_key": audit["join_key"],
    "n_manifest_samples": audit["n_manifest_samples"],
    "n_disagreements": audit["n_disagreements"],
    "manifest_sha256": audit["manifest_sha256"],
    "release_state": "candidate_only",
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Clinical provenance freeze: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
