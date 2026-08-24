#!/usr/bin/env bash
#SBATCH --job-name=masld-mpra-outcomes
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ANALYSIS="${PROJECT}/GWAS/finemapping/results/seqfunc/mpra_benchmark/v2"
HEPG2="${ANALYSIS}/source_reproduction/HepG2.allele_replicate_counts.tsv.gz"
HEPG2_SHA=9ccddf507a79b9cc3ccc3129479cea047fa74c32f66c45d91279a41287c68a26
LX2="${ANALYSIS}/source_reproduction/LX2.allele_replicate_counts.tsv.gz"
LX2_SHA=69f66f07247b3c2a0c018ac19476803b2896e697d89686cf4b4ea55a77dfb9a9
GATE="${ANALYSIS}/gate_verdict.json"
GATE_SHA=aeb08120ab28884cef3fcd1ff86985e8420c0fdebccf13ca4dd2e4b4f7c3586d
MAPPING="${PROJECT}/GWAS/finemapping/data/seqfunc_external/hu2025_mpra/raw_counts_v2/GSE281364_Lx2_Ctrl_pGL424_NAFLD-CRS-bc_pair_filtered_rep1_DNA_RNA_counts.txt.gz"
MAPPING_SHA=7d29c639cc9c59a59b9a74f55421cb2aabc617abe298f657050f688acbd1a57a
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
TARGET="${ROOT}/executions/gse281364-replicate-outcomes-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed MPRA outcome stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite MPRA outcome artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n%s  %s\n%s  %s\n%s  %s\n' \
  "${HEPG2_SHA}" "${HEPG2}" \
  "${LX2_SHA}" "${LX2}" \
  "${GATE_SHA}" "${GATE}" \
  "${MAPPING_SHA}" "${MAPPING}" \
  | sha256sum --check --strict

"${PYTHON}" "${ROOT}/scripts/audit_gse281364_replicate_outcomes.py" \
  --hepg2-counts "${HEPG2}" \
  --lx2-counts "${LX2}" \
  --mapping-counts "${MAPPING}" \
  --gate "${GATE}" \
  --output "${STAGE}/outcomes" \
  > "${STAGE}/validation/audit.stdout.txt" \
  2> "${STAGE}/validation/audit.stderr.txt"
"${PYTHON}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/audit_gse281364_replicate_outcomes.py" \
  "${ROOT}/slurm/gse281364_admit_replicate_outcomes.sbatch" \
  "${HEPG2}" "${LX2}" "${MAPPING}" "${GATE}" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "outcomes/outcome_admission_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["rows"] != 172_752
    or receipt["elements"] != 5_442
    or receipt["paired_construct_elements"] != 5_355
    or receipt["reference_only_construct_elements"] != 81
    or receipt["alternative_only_construct_elements"] != 6
    or receipt["observed_rows"] != 172_528
    or receipt["below_qc_rows"] != 224
    or receipt["missing_evidence_encoded_as_zero"]
    or receipt["experimental_replicates"] != 16
    or receipt["donor_count"] != 0
    or receipt["replicates_are_independent_donors"]
    or receipt["activity_transform_fit_scope"] != "not_fit_source_admission_only"
    or receipt["eQTL_or_ieQTL_supervision"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["champion_eligible"]
):
    raise SystemExit("MPRA replicate outcome receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_replicate_safe_outcomes",
        "hepg2_source_sha256": "9ccddf507a79b9cc3ccc3129479cea047fa74c32f66c45d91279a41287c68a26",
        "lx2_source_sha256": "69f66f07247b3c2a0c018ac19476803b2896e697d89686cf4b4ea55a77dfb9a9",
        "source_reproduction_gate_sha256": "aeb08120ab28884cef3fcd1ff86985e8420c0fdebccf13ca4dd2e4b4f7c3586d",
        "mapping_source_sha256": "7d29c639cc9c59a59b9a74f55421cb2aabc617abe298f657050f688acbd1a57a",
        "contexts": 4,
        "experimental_replicates": 16,
        "observed_rows": 172528,
        "explicit_below_qc_rows": 224,
        "missing_evidence_encoded_as_zero": False,
        "donor_count": 0,
        "outcome_role": "exposed_development_MPRA_only",
        "sealed_outcomes_loaded": False,
        "champion_eligible": False,
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 replicate-safe outcomes: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
