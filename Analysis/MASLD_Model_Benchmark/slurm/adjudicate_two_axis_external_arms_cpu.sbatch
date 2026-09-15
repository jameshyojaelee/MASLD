#!/usr/bin/env bash
#SBATCH --job-name=model-check-952
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

BULK=${PROJECT}/Analysis/Deconvolution/bulk
PRIOR=${ROOT}/config/evaluation/bulk_existing_source_pool_composition_activation.json
CONTRACT=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
INSTALLED=${ROOT}/config/evaluation/two_axis_activity_fibrosis_external_activation.json

TARGET=${ROOT}/executions/model-check-952-${SLURM_JOB_ID}-two-axis-activation
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite: %s\n' "${TARGET}" >&2; exit 1
fi
if [[ -e "${INSTALLED}" ]]; then
  printf 'Refusing to overwrite an existing activation record: %s\n' "${INSTALLED}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
hostname; "${PYTHON}" -V
bash -n slurm/adjudicate_two_axis_external_arms_cpu.sbatch

# Every input checked before a byte is read. The two frozen bases are read-only
# to this step and are re-checked afterwards to prove they were not edited.
cat > "${STAGE}/validation/inputs.sha256" <<'SHA'
6fb3cccd8c79d3e445725329086f8adfe1904753d2fd3c5c975be0ec06037de1  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/config/evaluation/bulk_existing_source_pool_composition_activation.json
b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/config/evaluation/microarray_transfer_label_semantics.toml
8e8e0661c03a431f7006cf06c1e9bfdf62e4030392c505a4fd93e5cda1a16cba  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE130970/GSE130970_metadata.tsv
0d9aa05ce52b7c16501a947e12e25a9b0513a83a2a3a26751307de979a500334  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE130970/GSE130970_counts.tsv
667782f577e8ad7cc3d80d5b776d02a4e32c20e8fe54a964ad3102069e1ed7d9  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE193066/GSE193066_metadata.tsv
a1b06183e9c9f4444db264afcf2fd18d0215a3d68a7163fd3bc0d3851acd3667  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE193066/GSE193066_counts.tsv
f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json
SHA
sha256sum --check --strict "${STAGE}/validation/inputs.sha256"

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_adjudicate_two_axis_external_arms \
  tests.unit.test_build_showcase_aspect_lineage_substrate \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" \
  scripts/adjudicate_two_axis_external_arms.py \
  --bulk-root "${BULK}" \
  --prior-activation "${PRIOR}" \
  --label-contract "${CONTRACT}" \
  --family-sizes 4 12 \
  --output "${STAGE}/activation" \
  > "${STAGE}/build_stdout.json" 2> "${STAGE}/build_stderr.txt"

sha256sum --check --strict "${STAGE}/validation/inputs.sha256"
printf 'frozen bases and cohort metadata byte-identical before and after\n'

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/adjudicate_two_axis_external_arms.py \
  scripts/build_showcase_aspect_lineage_substrate.py \
  tests/unit/test_adjudicate_two_axis_external_arms.py \
  slurm/adjudicate_two_axis_external_arms_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
act = root / "activation"
r = json.loads((act / "two_axis_activity_fibrosis_external_activation.json").read_text())
a = r["adjudications"]["GSE130970"]
b = r["adjudications"]["GSE193066"]
de = r["detectable_effect"]

if (
    r["frozen_bases_edited"] is not False
    or a["unit_census"]["metadata_rows"] != 78
    or a["unit_census"]["distinct_people_assertable"] is not False
    or a["fibrosis_standard_assessment"]["all_five_levels_distinguishable"] is not True
    or a["fibrosis_standard_assessment"]["verdict"] != "admissible_for_the_graded_primary"
    or a["activity_assessment"]["component_sum_mismatches"] != 0
    or a["both_axes_on_every_row"] is not True
    or a["control_field_adjudication"]["control_or_disease_field_present"] is not False
    or b["unit_census"]["metadata_rows"] != 164
    or b["unit_census"]["distinct_participants"] != 106
    or b["unit_census"]["orphan_second_biopsies"] != []
    or b["fibrosis_standard_assessment"]["all_five_levels_distinguishable"] is not False
    or b["fibrosis_standard_assessment"]["omitted_level"] != 4
    or "separate_omitted_level_arm" not in b["fibrosis_standard_assessment"]["verdict"]
    or r["separate_arm_rule"]["arms_are_never_pooled"] is not True
    or r["separate_arm_rule"]["arms_are_never_meta_analysed_against_each_other"] is not True
    or r["source_training_scope_ruling"]["ruling"] != "out_of_scope"
    or r["source_training_scope_ruling"]["contract_quotes_verified_present"] is not True
    or r["approval"]["technical_adjudication"]["user_approved"] is not False
    or r["approval"]["prior_activation_approval_does_not_extend"]["extends_to_this_record"] is not False
    or de["graded_primary_fibrosis"]["arm_can_detect_an_effect_below_its_ceiling"] is not True
    or de["omitted_level_arm_fibrosis"]["arm_can_detect_an_effect_below_its_ceiling"] is not True
):
    raise SystemExit("two-axis activation verification differs")

freeze_tree(act, {
    "artifact_class": "two_axis_activity_fibrosis_external_activation",
    "gse130970_verdict": a["fibrosis_standard_assessment"]["verdict"],
    "gse193066_verdict": b["fibrosis_standard_assessment"]["verdict"],
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(act)
freeze_tree(root, {
    "artifact_class": "two_axis_activity_fibrosis_external_activation_campaign",
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "GSE130970": {
        "n": a["unit_census"]["metadata_rows"],
        "verdict": a["fibrosis_standard_assessment"]["verdict"],
        "fibrosis_ceiling": de["graded_primary_fibrosis"]["attainable_ceiling_rho"],
        "fibrosis_mde": de["graded_primary_fibrosis"]["minimum_detectable_rho_at_80_percent_power"],
    },
    "GSE193066": {
        "participants": b["unit_census"]["distinct_participants"],
        "verdict": b["fibrosis_standard_assessment"]["verdict"],
        "fibrosis_ceiling": de["omitted_level_arm_fibrosis"]["attainable_ceiling_rho"],
        "fibrosis_mde": de["omitted_level_arm_fibrosis"]["minimum_detectable_rho_at_80_percent_power"],
    },
}, indent=2, sort_keys=True))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
# Install the record beside the prior activation, refusing to clobber.
cp -n "${TARGET}/activation/two_axis_activity_fibrosis_external_activation.json" "${INSTALLED}"
printf 'activation: %s\n' "${TARGET}"
printf 'installed record: %s\n' "${INSTALLED}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/activation/ARTIFACTS.json" "${INSTALLED}"
