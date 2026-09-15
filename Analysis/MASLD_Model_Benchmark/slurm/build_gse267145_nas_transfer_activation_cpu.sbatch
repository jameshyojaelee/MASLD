#!/usr/bin/env bash
#SBATCH --job-name=model-data-929
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-880-21130257-gse135251-source
SOURCE_SHA=4dab6a759371ccca56227bc3ce1df2c4e9ad27397cab07db4c337669c5003e8d
TARGET=${ROOT}/executions/model-data-064-21079902/fixture/molecular
TARGET_SHA=2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456
AUDIT=${ROOT}/executions/model-data-928-21132819/audit
AUDIT_SHA=2de5a671beb9a4d596700e7ba430908d9d86fa9f85df8989660a629acac78da6
TASK_SPEC=${ROOT}/config/evaluation/gse267145_nas_sum_transfer_task.toml
TASK_SPEC_SHA=5c1505d1b1fa05f3513823573146a8372cd0a8cd6e6c3e95605adb7ec09f3c78
GATE=${ROOT}/config/evaluation/gse267145_nas_sum_transfer_promotion_gate.json
GATE_SHA=09b900c4ad8167a0d2cb4fa5ce49fdd6a06fc918bd22cbba313deb154bd42c51
OUT=${ROOT}/executions/model-data-929-${SLURM_JOB_ID}
STAGE=${OUT}.staging
FAILED=${OUT}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${OUT}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite NAS transfer activation artifacts.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_SHA}" "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${TARGET_SHA}" "${TARGET}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${AUDIT_SHA}" "${AUDIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${TASK_SPEC_SHA}" "${TASK_SPEC}" | sha256sum --check --strict
printf '%s  %s\n' "${GATE_SHA}" "${GATE}" | sha256sum --check --strict
bash -n slurm/build_gse267145_nas_transfer_activation_cpu.sbatch
# The TaskSpec pins its own check, and the contamination audit must have passed
# before any model in this lane is allowed to exist.
"${PYTHON}" - "${TASK_SPEC}" "${GATE_SHA}" "${AUDIT}/audit.json" <<'PY'
import hashlib, json, re, sys
text = open(sys.argv[1], encoding="utf-8").read()
m = re.search(r'promotion_gate_config_sha256\s*=\s*"([0-9a-f]{64})"', text)
if m is None or m.group(1) != sys.argv[2]:
    raise SystemExit("TaskSpec does not pin this promotion gate")
audit = json.load(open(sys.argv[3], encoding="utf-8"))
if audit["status"] != "pass_no_overlap_detected":
    raise SystemExit(f"contamination audit did not pass: {audit['status']}")
if audit["identifier_layer"]["shared_identifiers_total"] != 0:
    raise SystemExit("contamination audit reports shared identifiers")
if audit["molecular_layer"]["reciprocal_best_hits_exceeding_reference"]:
    raise SystemExit("contamination audit reports a molecular near-duplicate")
print("TaskSpec gate pin and contamination audit both verified")
PY
"${PYTHON}" -m unittest -v \
  tests.unit.test_build_gse267145_nas_transfer_activation \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/build_gse267145_nas_transfer_activation.py \
  --source-feature-axis "${SOURCE}/source/molecular/rna_feature_axis.tsv" \
  --target-feature-axis "${TARGET}/rna_feature_axis.tsv" \
  --target-participant-axis "${TARGET}/participant_axis.tsv" \
  --expected-source-genes 73720 \
  --expected-target-genes 42163 \
  --expected-target-participants 99 \
  --output "${STAGE}/activation" \
  > "${STAGE}/activation_stdout.json" 2> "${STAGE}/activation_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/build_gse267145_nas_transfer_activation.py \
  tests/unit/test_build_gse267145_nas_transfer_activation.py \
  slurm/build_gse267145_nas_transfer_activation_cpu.sbatch \
  "${TASK_SPEC}" "${GATE}" "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" "${TARGET}/ARTIFACTS.json" "${AUDIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv, json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
act = root / "activation"
r = json.loads((act / "activation_receipt.json").read_text(encoding="utf-8"))
with (act / "common_stable_gene_axis.tsv").open(encoding="utf-8", newline="") as h:
    genes = [row["stable_gene_id"] for row in csv.DictReader(h, delimiter="\t")]
with (act / "target_participant_roster.tsv").open(encoding="utf-8", newline="") as h:
    roster = list(csv.DictReader(h, delimiter="\t"))
forbidden = {"nash_crn_component_sum", "steatosis", "ballooning",
             "lobular_inflammation", "lobular_necrosis", "fibrosis", "stage3", "stage5"}
if (
    r.get("status") != "pass_label_blind_shared_axis_and_target_roster"
    or r.get("target_outcome_columns_read") is not False
    or r.get("expression_values_read") is not False
    or r.get("shared_genes") != len(genes)
    or r.get("target_participants") != 99
    or genes != sorted(genes) or len(set(genes)) != len(genes)
    or len(roster) != 99
    or forbidden & set(roster[0])
):
    raise SystemExit("NAS transfer activation verification differs")
freeze_tree(act, {
    "artifact_class": "gse267145_nas_sum_transfer_activation",
    "shared_genes": len(genes), "target_participants": 99,
    "target_outcome_columns_read": False, "status": "passed",
})
verify_frozen_tree(act)
freeze_tree(root, {
    "artifact_class": "gse267145_nas_sum_transfer_activation_campaign",
    "target_outcomes_read": False, "metrics_calculated": False, "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({"shared_genes": len(genes)}, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${OUT}"
sha256sum "${OUT}/ARTIFACTS.json" "${OUT}/activation/ARTIFACTS.json"
