#!/usr/bin/env bash
#SBATCH --job-name=model-check-1051
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Corrects the recorded *reason* the participant bootstrap was dropped.
#
# Additive overlay over Stage 0b's result and Stage 0c's prespecification.
# Neither is retrofitted: the Arm B prespecification is held back, and the Arm B
# analysis script was hashed into a completed job's source manifest, so editing
# it would record a digest for content that job did not execute.
#
# The correction carries its own evidence. It recomputes the df comparison, the
# plain-correlation bootstrap bias that bounds the rule's scope, and Stage 0's
# C3 intervals, rather than quoting numbers from a report.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/record_bootstrap_df_correction.py
ARM_A_PRESPEC=${ROOT}/executions/model-check-1046-21181811-axis-count-prespec/prespecification
ARM_A_RESULT=${ROOT}/executions/model-check-1047-21181822-axis-count/axis_count
ARM_B_PRESPEC=${ROOT}/executions/model-check-1049-21182303-arm-b-axis-count-prespec/prespecification
ENDPOINTS=${ROOT}/executions/model-data-061-21079623/activation/participant_endpoints.tsv
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular

TARGET=${ROOT}/executions/model-check-1051-${SLURM_JOB_ID}-bootstrap-df-correction
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the correction record.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/record_bootstrap_df_correction_cpu.sbatch

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
e7397c3038dfbbdd1619048acfdd6cabc0e018a7af161175925130a7a6ffed3e  executions/model-check-1046-21181811-axis-count-prespec/prespecification/ARTIFACTS.json
78ce955b6ef43ff45ebf78ec7744aa464febf3da93ce8aff4b6daf0f80c4ea3e  executions/model-check-1047-21181822-axis-count/axis_count/ARTIFACTS.json
91b18fbe8cb061f04825beb3d9c7651ade92e6055782eacd06b6ece9c2cab49f  executions/model-check-1048-21182175-axis-count-derived-after/addendum/ARTIFACTS.json
58f72a4522309f9d4e0dfb1630c37fd994a308e714856b88508aa40cb27732f6  executions/model-check-1049-21182303-arm-b-axis-count-prespec/prespecification/ARTIFACTS.json
9a13c0fd65db0483566b7495def11ade3339861de091ea29fcbdc4c34f47b5b9  executions/model-data-061-21079623/activation/participant_endpoints.tsv
812d5fce6e7cc88f6939ae493db7b7e7470057c21cf3b219abc39f587ff64ddc  executions/model-data-064-21079902/fixture/molecular/rna_values.npy
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m unittest -v tests.unit.test_record_bootstrap_df_correction \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" \
  --arm-a-prespec "${ARM_A_PRESPEC}" \
  --arm-a-result "${ARM_A_RESULT}" \
  --arm-b-prespec "${ARM_B_PRESPEC}" \
  --endpoints "${ENDPOINTS}" \
  --molecular "${MOLECULAR}" \
  --output "${STAGE}/correction" \
  --genes 6000 \
  --resamples 400 \
  --label-resamples 10000 \
  --seed 20260827 \
  > "${STAGE}/correction_stdout.json" 2> "${STAGE}/correction_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" \
  slurm/record_bootstrap_df_correction_cpu.sbatch \
  tests/unit/test_record_bootstrap_df_correction.py \
  "${ENDPOINTS}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
payload = json.loads(
    (root / "correction" / "bootstrap_df_correction.json").read_text(encoding="utf-8")
)
if payload["record_is_additive_overlay"] is not True:
    raise SystemExit("the record does not declare itself an overlay")
for key in ("is_a_gate", "frozen_verdict_altered"):
    if payload[key] is not False:
        raise SystemExit(f"{key} must be false for a correction record")
if payload["stage_0b_outcome_unchanged"] != "ONE_AXIS_ACTIVITY":
    raise SystemExit("the Stage 0b outcome moved")
if payload["stage_0c_design_unchanged"] is not True:
    raise SystemExit("the Stage 0c design moved")
if payload["what_is_corrected"]["only_the_reason_was_wrong"] is not True:
    raise SystemExit("the record claims more than the reason was wrong")
if payload["honesty_about_ordering"]["may_not_be_cited_as_pre_registered"] is not True:
    raise SystemExit("the overlay does not forbid citing it as pre-registered")

evidence = payload["evidence_measured_here"]
# The correction's central claim must hold in its own measured evidence.
if not (
    evidence["bootstrap_median_count_at_df_96"]
    > evidence["bootstrap_median_count_at_df_distinct_minus_3"]
):
    raise SystemExit("the measured evidence does not support the correction")
if not 50.0 < evidence["mean_distinct_participants_per_resample"] < 75.0:
    raise SystemExit("the duplication rate is not in the expected range for n=99")

scope = payload["scope_of_the_rule"]
for entry in scope["plain_correlation_bootstraps_are_unaffected"]["per_pair"]:
    if abs(entry["bias_median_minus_observed"]) > 0.01:
        raise SystemExit(
            "a plain-correlation bootstrap shows material bias; the rule's "
            "stated scope is wrong"
        )
c3 = scope["stage_0_c3_association_vector_intervals"]
if c3["largest_downward_bias"] > 0.0:
    raise SystemExit(
        "the C3 bias is not downward; the reported reading of Stage 0 must change"
    )
if c3["no_stage_0_conclusion_moves"] is not True:
    raise SystemExit("a Stage 0 conclusion moved")

verify_frozen_tree(root / "correction")
freeze_tree(root, {
    "artifact_class": "masld_bench_bootstrap_df_correction_campaign",
    "record_is_additive_overlay": True,
    "corrects": "a recorded explanation, not a statistic",
    "frozen_verdict_altered": False,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "observed_count_at_df_96": evidence["observed_count_at_df_96"],
    "mean_distinct_participants": evidence["mean_distinct_participants_per_resample"],
    "bootstrap_median_at_df_96": evidence["bootstrap_median_count_at_df_96"],
    "bootstrap_median_at_honest_df": evidence[
        "bootstrap_median_count_at_df_distinct_minus_3"],
    "largest_c3_downward_bias": c3["largest_downward_bias"],
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/correction/ARTIFACTS.json"
