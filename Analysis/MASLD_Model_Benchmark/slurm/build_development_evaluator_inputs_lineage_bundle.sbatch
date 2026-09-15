#!/usr/bin/env bash
#SBATCH --job-name=seq-eval-inputs
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${LINEAGE:?LINEAGE is required}"
: "${OUTCOME_BUNDLE_JOB_ID:?OUTCOME_BUNDLE_JOB_ID is required}"
case "${LINEAGE}" in cholangiocyte|fibroblast|hepatocyte|macrophage|t_cell) ;;
  *) printf 'Invalid primary lineage: %s\n' "${LINEAGE}" >&2; exit 1 ;;
esac
[[ "${OUTCOME_BUNDLE_JOB_ID}" =~ ^[0-9]+$ ]] || { printf 'Invalid outcome bundle job ID\n' >&2; exit 1; }

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
DONOR_BIGWIGS="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
PSEUDOBULK="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
OUTCOME_BUNDLE="${ROOT}/executions/sequence-development-outcomes-${LINEAGE}-bundle-${OUTCOME_BUNDLE_JOB_ID}"
BUNDLE_TARGET="${ROOT}/executions/sequence-evaluator-inputs-${LINEAGE}-bundle-${SLURM_JOB_ID}"
BUNDLE_STAGE="${BUNDLE_TARGET}.staging"
BUNDLE_FAILED="${BUNDLE_TARGET}.failed"
preserve_bundle() { if [[ -d "${BUNDLE_STAGE}" && ! -e "${BUNDLE_FAILED}" ]]; then mv "${BUNDLE_STAGE}" "${BUNDLE_FAILED}"; fi; }
trap preserve_bundle EXIT
if [[ -e "${BUNDLE_TARGET}" || -e "${BUNDLE_STAGE}" || -e "${BUNDLE_FAILED}" ]]; then
  printf 'Refusing to overwrite evaluator-input bundle: %s\n' "${BUNDLE_TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${BUNDLE_STAGE}/validation" "${BUNDLE_STAGE}/receipts"
for pair in \
  "b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 ${BASE}/ARTIFACTS.json" \
  "96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa ${IMAGE}" \
  "61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 ${CORE}/ARTIFACTS.json" \
  "01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f ${DONOR_BIGWIGS}/ARTIFACTS.json" \
  "3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b ${PSEUDOBULK}/ARTIFACTS.json" \
  "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f ${SPLIT}/ARTIFACTS.json" \
  "afbefe391bfbd3651569b2372f42fbe6a4bd30288d596e62883075f4bd78ba29 ${ROOT}/scripts/build_chrombpnet_development_profiles.py" \
  "a64296babe6bd1219f7832d42bc9cda005d49ceac7e56183cf3b84570a80088a ${ROOT}/tests/unit/test_build_chrombpnet_development_profiles.py" \
  "5f6c9125438eb8ca2ed1233731b4a34a4d717d2357c49fa905e8e06d96c21746 ${ROOT}/scripts/build_sequence_training_mean_baseline.py" \
  "960086c44001a27a3a1adb748d8ad9f4ed380041ebf87a6c26faf6250fce7297 ${ROOT}/tests/unit/test_build_sequence_training_mean_baseline.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
printf '%s  %s\n' 01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f \
  "${DONOR_BIGWIGS}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${OUTCOME_BUNDLE}" "${SPLIT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v \
    tests.unit.test_build_chrombpnet_development_profiles \
    tests.unit.test_build_sequence_training_mean_baseline \
  > "${BUNDLE_STAGE}/validation/tests.stdout.txt" \
  2> "${BUNDLE_STAGE}/validation/tests.stderr.txt"

build_fold() {
  local fold="$1" split_id outcome_receipt outcome_root outcome_sha
  local profile_target profile_stage baseline_target baseline_stage receipt
  split_id="donor${fold}_genomic${fold}"
  receipt="${BUNDLE_STAGE}/receipts/${split_id}.json"
  if [[ "${LINEAGE}" == hepatocyte && "${fold}" == 0 ]]; then
    local existing_profile existing_profile_sha existing_baseline existing_baseline_sha
    existing_profile="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-profiles-21064545"
    existing_profile_sha=a0cf1a5d62f4d8aaa78af8f5368e2ee7c33de99cea7867eeda90a74c7e335dba
    existing_baseline="${ROOT}/executions/sequence-training-mean-hepatocyte-donor0-genomic0-21064587"
    existing_baseline_sha=76980c7ba990aa5948feedcfd02f13ef3d3f550f702512b7a7fd1ccceca4fc82
    printf '%s  %s\n' "${existing_profile_sha}" "${existing_profile}/ARTIFACTS.json" | sha256sum --check --strict || return
    printf '%s  %s\n' "${existing_baseline_sha}" "${existing_baseline}/ARTIFACTS.json" | sha256sum --check --strict || return
    module purge && module load python/3.11.5-GCCcore-13.2.0 || return
    PYTHONPATH="${ROOT}/src" python - "${existing_profile}" "${existing_baseline}" \
      "${receipt}" "${split_id}" "${existing_profile_sha}" "${existing_baseline_sha}" <<'PY' || return
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree
profile, baseline, receipt = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
for root in (profile, baseline): verify_frozen_tree(root)
receipt.write_text(json.dumps({
    "schema_version": "masld-bench-sequence-evaluator-input-bundle-task-v1",
    "status": "reused_frozen_artifacts", "split_id": sys.argv[4],
    "profile_root": str(profile.resolve()), "profile_artifacts_sha256": sys.argv[5],
    "baseline_root": str(baseline.resolve()), "baseline_artifacts_sha256": sys.argv[6],
}, sort_keys=True, indent=2) + "\n")
PY
    return 0
  fi
  outcome_receipt="${OUTCOME_BUNDLE}/receipts/${split_id}.json"
  readarray -t outcome_values < <(python3 - "${outcome_receipt}" <<'PY'
import json,sys
x=json.load(open(sys.argv[1])); print(x["artifact_root"]); print(x["artifacts_sha256"])
PY
  )
  [[ "${#outcome_values[@]}" -eq 2 ]] || return
  outcome_root="${outcome_values[0]}"
  outcome_sha="${outcome_values[1]}"
  [[ "${outcome_root}" == "${ROOT}"/executions/* ]] || return
  printf '%s  %s\n' "${outcome_sha}" "${outcome_root}/ARTIFACTS.json" | sha256sum --check --strict || return
  profile_target="${ROOT}/executions/chrombpnet-${LINEAGE}-${split_id}-profiles-${SLURM_JOB_ID}"
  profile_stage="${profile_target}.staging"
  baseline_target="${ROOT}/executions/sequence-training-mean-${LINEAGE}-${split_id}-${SLURM_JOB_ID}"
  baseline_stage="${baseline_target}.staging"
  if [[ -e "${profile_target}" || -e "${profile_stage}" || -e "${baseline_target}" || -e "${baseline_stage}" ]]; then
    return 1
  fi
  install -d -m 0750 "${profile_stage}/validation" "${baseline_stage}/validation" || return
  module purge && module load singularity/4.2.2 || return
  singularity exec --cleanenv --containall --net --network none \
    --bind "${ROOT}:${ROOT}:rw" \
    --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
    --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
    "${IMAGE}" python "${ROOT}/scripts/build_chrombpnet_development_profiles.py" \
      --donor-bigwigs "${DONOR_BIGWIGS}" \
      --donor-bigwigs-artifacts-sha256 01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f \
      --split-contract "${SPLIT}" \
      --split-artifacts-sha256 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
      --regional-outcomes "${outcome_root}" \
      --regional-outcomes-artifacts-sha256 "${outcome_sha}" \
      --split-id "${split_id}" --lineage "${LINEAGE}" \
      --output "${profile_stage}/profiles" --workers 8 \
      > "${profile_stage}/validation/build.stdout.txt" \
      2> "${profile_stage}/validation/build.stderr.txt" || return
  singularity exec --cleanenv --containall --net --network none \
    --bind "${ROOT}:${ROOT}:rw" \
    --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
    --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
    "${IMAGE}" python "${ROOT}/scripts/build_sequence_training_mean_baseline.py" \
      --pseudobulk "${PSEUDOBULK}" \
      --pseudobulk-artifacts-sha256 3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b \
      --split-contract "${SPLIT}" \
      --split-artifacts-sha256 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
      --split-id "${split_id}" --lineage "${LINEAGE}" \
      --output "${baseline_stage}/baseline" \
      > "${baseline_stage}/validation/build.stdout.txt" \
      2> "${baseline_stage}/validation/build.stderr.txt" || return
  singularity exec --cleanenv --containall --net --network none \
    --bind "${ROOT}:${ROOT}:rw" --env PYTHONNOUSERSITE=1 \
    "${IMAGE}" python - "${profile_stage}" "${baseline_stage}" "${split_id}" "${LINEAGE}" <<'PY' || return
from pathlib import Path
import json
import sys
import h5py

profile, baseline, split_id, lineage = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
ps = json.loads((profile / "profiles/summary.json").read_text())
bs = json.loads((baseline / "baseline/summary.json").read_text())
if ps["split_id"] != split_id or ps["lineage_id"] != lineage or bs["split_id"] != split_id or bs["lineage_id"] != lineage:
    raise SystemExit("evaluator input identity differs")
if ps["model_inference_input_eligible"] or ps["benchmark_metrics_calculated"] or bs["held_donor_outcomes_exposed"] or bs["benchmark_metrics_calculated"]:
    raise SystemExit("evaluator input firewall differs")
if ps["roles"]["valid"]["windows"] != 16000 or ps["roles"]["test"]["windows"] != 16000 or bs["roles"] != {"valid": 16000, "test": 16000}:
    raise SystemExit("evaluator input window census differs")
with h5py.File(profile / "profiles/profiles.h5", "r") as handle:
    if handle["valid/counts"].shape[1:] != (16000, 1000): raise SystemExit("profile shape differs")
with h5py.File(baseline / "baseline/training_mean_baseline.h5", "r") as handle:
    if handle["valid/pseudobulk_counts"].shape != (16000, 1000): raise SystemExit("baseline shape differs")
PY
  sha256sum "${ROOT}/scripts/build_chrombpnet_development_profiles.py" \
    "${ROOT}/tests/unit/test_build_chrombpnet_development_profiles.py" \
    "${ROOT}/slurm/build_development_evaluator_inputs_lineage_bundle.sbatch" \
    > "${profile_stage}/source.sha256" || return
  sha256sum "${ROOT}/scripts/build_sequence_training_mean_baseline.py" \
    "${ROOT}/tests/unit/test_build_sequence_training_mean_baseline.py" \
    "${ROOT}/slurm/build_development_evaluator_inputs_lineage_bundle.sbatch" \
    > "${baseline_stage}/source.sha256" || return
  module purge && module load python/3.11.5-GCCcore-13.2.0 || return
  PYTHONPATH="${ROOT}/src" python - "${profile_stage}" "${baseline_stage}" \
    "${split_id}" "${LINEAGE}" "${outcome_sha}" "${receipt}" <<'PY' || return
from pathlib import Path
from hashlib import sha256
import json
import sys
from masld_bench.artifacts import freeze_tree

profile, baseline, split_id, lineage, outcome_sha, receipt = (
    Path(sys.argv[1]).resolve(strict=True), Path(sys.argv[2]).resolve(strict=True),
    sys.argv[3], sys.argv[4], sys.argv[5], Path(sys.argv[6])
)
freeze_tree(profile, {
    "artifact_class": "chrombpnet_development_base_profiles", "dataset_id": "gse296875",
    "split_id": split_id, "lineage_id": lineage, "biological_unit": "donor",
    "regional_outcomes_artifacts_sha256": outcome_sha, "observed_atac_outcomes": True,
    "model_inference_input_eligible": False, "benchmark_metrics_calculated": False,
    "regional_outcomes_reproduced": True, "status": "passed",
})
freeze_tree(baseline, {
    "artifact_class": "sequence_training_mean_baseline", "model_id": "training_pseudobulk_mean",
    "dataset_id": "gse296875", "split_id": split_id, "lineage_id": lineage,
    "biological_unit": "donor", "held_donor_outcomes_exposed": False,
    "benchmark_metrics_calculated": False, "status": "passed",
})
profile_sha = sha256((profile / "ARTIFACTS.json").read_bytes()).hexdigest()
baseline_sha = sha256((baseline / "ARTIFACTS.json").read_bytes()).hexdigest()
receipt.write_text(json.dumps({
    "schema_version": "masld-bench-sequence-evaluator-input-bundle-task-v1",
    "status": "passed", "split_id": split_id,
    "profile_root": str(profile).removesuffix(".staging"),
    "profile_artifacts_sha256": profile_sha,
    "baseline_root": str(baseline).removesuffix(".staging"),
    "baseline_artifacts_sha256": baseline_sha,
}, sort_keys=True, indent=2) + "\n")
PY
  mv "${profile_stage}" "${profile_target}" || return
  mv "${baseline_stage}" "${baseline_target}" || return
}

failures=0
for fold in 0 1 2 3 4; do
  if ! build_fold "${fold}"; then
    split_id="donor${fold}_genomic${fold}"
    for stage in \
      "${ROOT}/executions/chrombpnet-${LINEAGE}-${split_id}-profiles-${SLURM_JOB_ID}.staging" \
      "${ROOT}/executions/sequence-training-mean-${LINEAGE}-${split_id}-${SLURM_JOB_ID}.staging"; do
      [[ ! -d "${stage}" || -e "${stage%.staging}.failed" ]] || mv "${stage}" "${stage%.staging}.failed"
    done
    printf '{"schema_version":"masld-bench-sequence-evaluator-input-bundle-task-v1","status":"failed","split_id":"%s"}\n' \
      "${split_id}" > "${BUNDLE_STAGE}/receipts/${split_id}.json"
    failures=$((failures + 1))
  fi
done

module purge
module load singularity/4.2.2
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${BUNDLE_STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_chrombpnet_development_profiles.py" \
  "${ROOT}/scripts/build_sequence_training_mean_baseline.py" \
  "${ROOT}/slurm/build_development_evaluator_inputs_lineage_bundle.sbatch" \
  > "${BUNDLE_STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BUNDLE_STAGE}" "${LINEAGE}" "${failures}" "${OUTCOME_BUNDLE_JOB_ID}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root, lineage, failures, source_job = Path(sys.argv[1]).resolve(strict=True), sys.argv[2], int(sys.argv[3]), sys.argv[4]
receipts = [json.loads(path.read_text()) for path in sorted((root / "receipts").glob("*.json"))]
if len(receipts) != 5: raise SystemExit("evaluator-input receipt count differs")
freeze_tree(root, {
    "artifact_class": "sequence_evaluator_inputs_lineage_bundle", "dataset_id": "gse296875",
    "lineage_id": lineage, "logical_tasks": 5, "failures": failures,
    "source_outcome_bundle_job_id": source_job, "biological_unit": "donor",
    "observed_atac_outcomes": True, "model_inference_input_eligible": False,
    "benchmark_metrics_calculated": False,
    "status": "passed" if failures == 0 else "partial_failure",
})
PY
trap - EXIT
mv "${BUNDLE_STAGE}" "${BUNDLE_TARGET}"
sha256sum "${BUNDLE_TARGET}/ARTIFACTS.json"
if [[ "${failures}" -ne 0 ]]; then exit 2; fi
