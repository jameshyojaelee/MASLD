#!/usr/bin/env bash
#SBATCH --job-name=masld-mpra-validate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RAW="${ROOT}/executions/gse281364-raw-reconstruction-21066288"
RAW_SHA=d4c41023575f30e052dc3d1ed55c4e957c764d08802365a8bcc28c71a2d3660c
OUTCOMES="${ROOT}/executions/gse281364-replicate-outcomes-21066307"
OUTCOMES_SHA=cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db
REFERENCE=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz
REFERENCE_SHA=9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215
PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
TARGET="${ROOT}/executions/gse281364-validated-reconstruction-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
WORK="${TMPDIR:-/tmp}/masld-mpra-validate-${SLURM_JOB_ID}"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed MPRA validation stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" || -e "${WORK}" ]]; then
  printf 'Refusing to overwrite MPRA validation output\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${WORK}"
for binding in "${RAW_SHA} ${RAW}" "${OUTCOMES_SHA} ${OUTCOMES}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
printf '%s  %s\n' "${REFERENCE_SHA}" "${REFERENCE}" \
  | sha256sum --check --strict
gzip -t "${REFERENCE}"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RAW}" "${OUTCOMES}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("MPRA validation frozen inputs verified")
PY
module purge

gzip -cd "${REFERENCE}" > "${WORK}/GRCh38.p14.genome.fa"
"${PYTHON}" "${ROOT}/scripts/validate_gse281364_reconstruction.py" \
  --constructs "${RAW}/reconstruction/constructs.tsv" \
  --reconstruction-stats "${RAW}/reconstruction/reconstruction_stats.json" \
  --reference-fasta "${WORK}/GRCh38.p14.genome.fa" \
  --replicate-outcomes "${OUTCOMES}/outcomes/replicate_outcomes.tsv.gz" \
  --output "${STAGE}/validated" \
  > "${STAGE}/validation/validator.stdout.txt" \
  2> "${STAGE}/validation/validator.stderr.txt"
"${PYTHON}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/validate_gse281364_reconstruction.py" \
  "${ROOT}/slurm/gse281364_validate_reconstruction.sbatch" \
  "${RAW}/ARTIFACTS.json" "${OUTCOMES}/ARTIFACTS.json" "${REFERENCE}" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validated/reconstruction_receipt.json").read_text())
qualified = receipt["qualified_paired_snv_elements"]
below = receipt["paired_construct_elements_below_sequence_qc"]
if (
    receipt["status"] != "pass"
    or receipt["constructs"] != 10_797
    or receipt["paired_construct_elements"] != 5_355
    or qualified <= 0
    or qualified + below != 5_355
    or receipt["reference_only_construct_elements"] != 81
    or receipt["alternative_only_construct_elements"] != 6
    or receipt["raw_reconstruction_original_threshold_passing_constructs"] != 0
    or receipt["raw_reconstruction_original_status_used_for_admission"]
    or receipt["replicate_outcome_rows"] != 172_752
    or receipt["observed_replicate_outcome_rows"] != 172_528
    or receipt["below_qc_replicate_outcome_rows"] != 224
    or receipt["missing_evidence_encoded_as_zero"]
    or receipt["experimental_replicates"] != 16
    or receipt["replicates_are_independent_donors"]
    or receipt["donor_count"] != 0
    or receipt["eQTL_or_ieQTL_supervision"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["champion_eligible"]
):
    raise SystemExit(f"MPRA validation receipt differs: {receipt!r}")
freeze_tree(root, {
    "artifact_class": "gse281364_validated_mpra_reconstruction",
    "raw_reconstruction_artifacts_sha256": "d4c41023575f30e052dc3d1ed55c4e957c764d08802365a8bcc28c71a2d3660c",
    "replicate_outcomes_artifacts_sha256": "cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db",
    "reference_fasta_sha256": "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215",
    "qualified_paired_snv_elements": qualified,
    "paired_construct_elements_below_sequence_qc": below,
    "missing_evidence_encoded_as_zero": False,
    "donor_count": 0,
    "outcome_role": "exposed_development_MPRA_only",
    "sealed_outcomes_loaded": False,
    "champion_eligible": False,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 validated reconstruction: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
