#!/usr/bin/env bash
#SBATCH --job-name=data-stage-260
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# One GSE83452 summarization shard. SHARD_INDEX is passed with --export.
#
# Each array is summarized on its own, so shards are independent by
# construction; that is what model-check-091 established on these CEL files.
# Sentinel arrays appear in every shard so the merge can confirm that
# independent jobs on different nodes agree bitwise.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

: "${SHARD_INDEX:?SHARD_INDEX must be exported}"
: "${PLAN_DIR:?PLAN_DIR must be exported}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/microarray-r-bioc-21104442
RUNTIME_SHA256=8115e0966804dfc4b7bfac3f95016896d099d82fc51f3428fab1ec34182a7d07
FIREWALL=${ROOT}/executions/model-check-091-21109479-firewall/firewall
FIREWALL_SHA256=414b3d39bde46150fdf2b8205e29b6e85327ebb971f165551b95290455be055b
CROSSWALK=${ROOT}/executions/model-data-087-21111875-crosswalk/crosswalk
CROSSWALK_SHA256=c6bfe6c3d44b98fd612cf81ba48014fb32968a4d690b9424369c7bc23018048d
CROSSWALK_TSV=${CROSSWALK}/gpl16686_gencode_v49_crosswalk.tsv
RAW_AUDIT=${ROOT}/executions/model-data-065-21079931
RAW_AUDIT_SHA256=8b30a2265c421f2917ed886d99cfaa08ee8d874abd40e68c8c72d48c726781d1
CEL_TAR=${RAW_AUDIT}/raw/GSE83452_RAW.tar
CEL_MANIFEST=${RAW_AUDIT}/GSE83452_cel_manifest.tsv
CEL_MANIFEST_SHA256=a7c20245ed131b27c56b256c450e9337d5f32d51a54e1b12e14f5e78339988f1

SHARD_FILE=${PLAN_DIR}/shard_${SHARD_INDEX}.txt
TARGET=${ROOT}/executions/data-stage-260-${SLURM_JOB_ID}-shard${SHARD_INDEX}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
CEL_WORK=${TMPDIR:-/tmp}/masld_gse83452_${SLURM_JOB_ID}

cleanup() {
  rm -rf "${CEL_WORK}"
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap cleanup EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then printf 'Refusing to overwrite: %s\n' "${path}" >&2; exit 1; fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation" "${STAGE}/shard"
install -d -m 0700 "${CEL_WORK}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}"      "${RUNTIME}/ARTIFACTS.json" \
  "${FIREWALL_SHA256}"     "${FIREWALL}/ARTIFACTS.json" \
  "${CROSSWALK_SHA256}"    "${CROSSWALK}/ARTIFACTS.json" \
  "${RAW_AUDIT_SHA256}"    "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST_SHA256}" "${CEL_MANIFEST}" \
  | sha256sum --check --strict
if [[ ! -f "${SHARD_FILE}" ]]; then printf 'Shard file absent: %s\n' "${SHARD_FILE}" >&2; exit 1; fi

cd "${ROOT}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"

# Refuse to summarize unless the GPL16686 separation itself passed.
python3 - "${FIREWALL}" <<'PY'
from pathlib import Path
import json, sys
v = json.loads((Path(sys.argv[1]) / "firewall_verdict.json").read_text(encoding="utf-8"))
if (
    v["status"] != "pass_inductive_firewall_on_real_CEL"
    or v["platform_id"] != "GPL16686"
    or v["arrays_per_summarization_call"] != 1
    or v["held_array_summary_invariant_to_other_arrays_in_the_run"] is not True
    or v["perturbation_control_detected_for_every_target"] is not True
    or v["gpl570_evidence_carried_over"] is not False
):
    raise SystemExit("refusing to summarize without a proven GPL16686 firewall")
PY

# Stage only this shard's arrays.
mapfile -t SHARD_ACCESSIONS < "${SHARD_FILE}"
STAGE_ARGS=()
for accession in "${SHARD_ACCESSIONS[@]}"; do
  [[ -n "${accession}" ]] && STAGE_ARGS+=(--accession "${accession}")
done
python3 scripts/stage_microarray_cel_arrays.py \
  --tar "${CEL_TAR}" --manifest "${CEL_MANIFEST}" --output "${CEL_WORK}" \
  --series GSE83452 --expected-records 231 --platform-id GPL16686 \
  "${STAGE_ARGS[@]}" \
  --receipt "${STAGE}/shard/cel_staging_receipt.json" \
  > "${STAGE}/validation/cel_staging_stdout.json" 2> "${STAGE}/validation/cel_staging_stderr.txt"

module purge
module load R/4.5.2-mba
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/summarize_gse83452_gpl16686_shard.R" \
  "${CEL_WORK}" "${SHARD_FILE}" "${CROSSWALK_TSV}" "${SHARD_INDEX}" "${STAGE}/shard" \
  > "${STAGE}/validation/shard_stdout.txt" 2> "${STAGE}/validation/shard_stderr.txt"

rm -rf "${CEL_WORK}"
cp "${SHARD_FILE}" "${STAGE}/shard/shard_accessions.txt"

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
sha256sum \
  scripts/summarize_gse83452_gpl16686_shard.R \
  scripts/stage_microarray_cel_arrays.py \
  slurm/summarize_gse83452_gpl16686_shard_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" "${FIREWALL}/ARTIFACTS.json" \
  "${CROSSWALK}/ARTIFACTS.json" "${CEL_MANIFEST}" \
  > "${STAGE}/source.sha256"

if find "${STAGE}" -type f -name '*.CEL' -print -quit | grep -q .; then
  printf 'Shard stage contains a raw CEL file\n' >&2; exit 1
fi

python3 - "${STAGE}" "${SHARD_INDEX}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
stage = Path(sys.argv[1]).resolve(strict=True)
shard_index = sys.argv[2]
shard = stage / "shard"
receipt = json.loads((shard / f"shard_{shard_index}_receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_label_blind_single_array_shard"
    or receipt["arrays_read_per_summarization_call"] != 1
    or receipt["matrices_round_trip_bitwise"] is not True
    or receipt["labels_read"] is not False
    or receipt["model_training_activated"] is not False
):
    raise SystemExit("shard receipt differs from the frozen contract")
metadata = {
    "artifact_class": "gse83452_gpl16686_single_array_shard",
    "series": "GSE83452", "platform_id": "GPL16686",
    "cohort_family_id": "antwerp_inserm_shared",
    "shard": shard_index,
    "arrays_summarized": receipt["arrays_summarized"],
    "arrays_read_per_summarization_call": 1,
    "matrices_round_trip_bitwise": True,
    "labels_read": False, "model_training_activated": False,
    "sealed_outcomes_read": False, "status": "passed",
}
freeze_tree(shard, metadata)
verify_frozen_tree(shard)
freeze_tree(stage, metadata | {"artifact_class": "gse83452_gpl16686_shard_campaign"})
verify_frozen_tree(stage)
PY

trap - EXIT
rm -rf "${CEL_WORK}"
mv "${STAGE}" "${TARGET}"
printf 'GSE83452 shard %s: %s\n' "${SHARD_INDEX}" "${TARGET}"
