#!/usr/bin/env bash
#SBATCH --job-name=model-check-898
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
TASKREC=${ROOT}/drafts/w3_prioritization/model_output_structural_corroboration_task.toml
TASKREC_SHA256=3b6284b354c24c35b85afd8f6096d57213c0f683316b1afea4c32fee25df5fcc
TEG=${ROOT}/config/tasks/typed_evidence_graph.toml
TEG_SHA256=c1237c32941a2b51c13a8e20741bb6697e541115de6c2c607964724c948bbb92
TARGET=${ROOT}/executions/model-check-898-${SLURM_JOB_ID}-w3-recorded-negative
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the W3 recorded negative.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/freeze_w3_recorded_negative_cpu.sbatch

# typed_evidence_graph.toml is read-only to this lane. Recorded BEFORE.
sha256sum "${TEG}" > "${STAGE}/typed_evidence_graph_before.sha256"
printf '%s  %s\n' "${TEG_SHA256}" "${TEG}" | sha256sum --check --strict
printf '%s  %s\n' "${TASKREC_SHA256}" "${TASKREC}" | sha256sum --check --strict

install -d -m 0750 "${STAGE}/record"
cp -p "${TASKREC}" "${STAGE}/record/w3_recorded_negative.toml"

# No check is registered for this lane, by ruling. Assert none was copied.
if ls "${STAGE}/record" | grep -qi 'gate'; then
  printf 'A gate file reached the frozen record; this lane must have none.\n' >&2; exit 1
fi

# AFTER.
sha256sum "${TEG}" > "${STAGE}/typed_evidence_graph_after.sha256"
if ! cmp -s "${STAGE}/typed_evidence_graph_before.sha256" "${STAGE}/typed_evidence_graph_after.sha256"; then
  printf 'typed_evidence_graph.toml changed during the freeze.\n' >&2; exit 1
fi
printf '%s  %s\n' "${TEG_SHA256}" "${TEG}" | sha256sum --check --strict

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum slurm/freeze_w3_recorded_negative_cpu.sbatch "${TASKREC}" "${TEG}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys, tomllib, json
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
rec = root / "record"
d = tomllib.loads((rec / "w3_recorded_negative.toml").read_bytes().decode())
if (d["status"] != "draft_substrate_cannot_support_this_lane_as_scoped"
        or d["verdict"]["finding"] != "The substrate cannot support this lane as scoped."
        or d["verdict"]["hepg2_only_realised"]["distinct_elements"] != 57
        or d["verdict"]["hepg2_only_realised"]["dav_distinct_rsid"] != 5
        or d["cell_type_specificity_evidence"]["dav_in_both"] != 13
        or d["lessons"]["convenience_is_a_tell"]["verbatim"] != "The convenience of doubling n was exactly the tell."
        or d["lessons"]["row_count_is_not_a_unit_count"]["first_class_finding"] is not True
        or d["lanes"]["region"]["status"] != "NO_ELIGIBLE_MODEL_OUTPUT"
        or "no_gate_registered" not in d
        or "campaign_level_finding" not in d):
    raise SystemExit("W3 recorded-negative verification differs")
freeze_tree(rec, {
    "artifact_class": "w3_prioritization_recorded_negative",
    "finding": "the substrate cannot support this lane as scoped",
    "promotion_gate_registered": False,
    "typed_evidence_graph_unmodified": True,
    "external_development_only": True,
    "champion_eligible": False,
    "status": "passed",
})
verify_frozen_tree(rec)
freeze_tree(root, {"artifact_class": "w3_recorded_negative_campaign", "status": "passed"})
verify_frozen_tree(root)
print(json.dumps({"finding": d["verdict"]["finding"],
                  "elements": d["verdict"]["hepg2_only_realised"]["distinct_elements"],
                  "dav_rsids": d["verdict"]["hepg2_only_realised"]["dav_distinct_rsid"]}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/record/ARTIFACTS.json"
cat "${TARGET}/typed_evidence_graph_after.sha256"
