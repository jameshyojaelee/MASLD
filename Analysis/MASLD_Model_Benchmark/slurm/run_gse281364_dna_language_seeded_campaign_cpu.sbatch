#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-train-605
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=36:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_dna_language_seeded_head_campaign.json
PROJECTION_CONFIG=${ROOT}/config/gse281364_sequence_context_refit_dag.toml
ROW_UNIVERSE=${ROOT}/executions/model-check-220-21088696/contract/row_universe.tsv
TARGET=${ROOT}/executions/model-cpu-train-605-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite DNA-language seeded campaign: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src
cd "${ROOT}"

SOURCE_FILES=(
  config/gse281364_dna_language_seeded_head_campaign.json
  config/gse281364_sequence_context_refit_dag.toml
  scripts/project_gse281364_long_range_embeddings.py
  scripts/gse281364_dna_lm_native_contract.py
  scripts/fit_gse281364_enformer_sei_heads.py
  scripts/fit_gse281364_open_sequence_heads.py
  scripts/evaluate_gse281364_dna_lm_common_lane.py
  scripts/reconcile_gse281364_dna_language_seeded_features.py
  scripts/fit_gse281364_dna_language_seeded_heads.py
  scripts/verify_gse281364_shared_control_identity.py
  scripts/evaluate_gse281364_dna_language_seeded_heads.py
  tests/unit/test_gse281364_dna_language_seeded_campaign.py
  slurm/run_gse281364_dna_language_seeded_campaign_cpu.sbatch
)
sha256sum "${SOURCE_FILES[@]}" > "${STAGE}/source.start.sha256"

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse281364_dna_language_seeded_campaign \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

# Complete the only missing long-range projection without opening outcomes.
"${PYTHON}" "${ROOT}/scripts/project_gse281364_long_range_embeddings.py" \
  --root "${ROOT}" \
  --config "${PROJECTION_CONFIG}" \
  --row-universe "${ROW_UNIVERSE}" \
  --model-roster nucleotide_transformer \
  --output "${STAGE}/nt_projection" \
  > "${STAGE}/validation/nt_projection.stdout.json" \
  2> "${STAGE}/validation/nt_projection.stderr.txt"

"${PYTHON}" - "${STAGE}/nt_projection" <<'PY'
import json
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_outcome_blind_long_range_projection"
    or receipt["models"] != ["nucleotide_transformer"]
    or receipt["elements"] != 1033
    or receipt["long_range_blocks"] != 239
    or receipt["outcomes_read"]
    or receipt["prediction_values_read"]
    or receipt["downstream_head_fit"]
    or receipt["sealed_assets_read"]
):
    raise SystemExit("Nucleotide Transformer projection firewall differs")
freeze_tree(root, {
    "artifact_class": "gse281364_nucleotide_transformer_long_range_projection",
    "dataset_id": "gse281364",
    "model_id": "nucleotide_transformer",
    "elements": 1033,
    "long_range_blocks": 239,
    "outcomes_read": False,
    "model_fit": False,
    "sealed_assets_read": False,
    "status": "pass_outcome_blind_long_range_projection",
})
verify_frozen_tree(root)
PY
NT_SHA=$(sha256sum "${STAGE}/nt_projection/ARTIFACTS.json" | awk '{print $1}')

# Freeze exact feature, checkpoint, runtime, exposure, missing piece, and TaskSpec status.
"${PYTHON}" "${ROOT}/scripts/reconcile_gse281364_dna_language_seeded_features.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --nt-projection-tree "${STAGE}/nt_projection" \
  --nt-projection-sha256 "${NT_SHA}" \
  --output "${STAGE}/task_spec" \
  > "${STAGE}/validation/reconciliation.stdout.json" \
  2> "${STAGE}/validation/reconciliation.stderr.txt"

"${PYTHON}" - "${STAGE}/task_spec" "${NT_SHA}" <<'PY'
import json
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_complete_feature_reconciliation_and_prespecified_taskspec"
    or receipt["feature_complete_models"] != ["caduceus", "dnabert2", "hyenadna", "nucleotide_transformer"]
    or receipt["feature_incomplete_models"]
    or receipt["features_imputed"]
    or receipt["obsolete_6kb_folds_reused"]
    or receipt["nucleotide_transformer_projection_artifacts_sha256"] != sys.argv[2]
    or receipt["outcomes_read"]
    or receipt["model_fit"]
    or receipt["sealed_assets_read"]
):
    raise SystemExit("DNA-language TaskSpec firewall differs")
freeze_tree(root, {
    "artifact_class": "gse281364_dna_language_feature_reconciliation_and_taskspec",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "models": ["caduceus", "dnabert2", "hyenadna", "nucleotide_transformer"],
    "candidate_count": 11,
    "elements": 1033,
    "long_range_blocks": 239,
    "feature_incomplete_models": [],
    "features_imputed": False,
    "obsolete_6kb_folds_reused": False,
    "outcomes_read": False,
    "model_fit": False,
    "sealed_assets_read": False,
    "status": "pass_complete_feature_reconciliation_and_prespecified_taskspec",
})
verify_frozen_tree(root)
PY
TASKSPEC_SHA=$(sha256sum "${STAGE}/task_spec/ARTIFACTS.json" | awk '{print $1}')

"${PYTHON}" "${ROOT}/scripts/fit_gse281364_dna_language_seeded_heads.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --task-spec-tree "${STAGE}/task_spec" \
  --task-spec-sha256 "${TASKSPEC_SHA}" \
  --nt-projection-tree "${STAGE}/nt_projection" \
  --nt-projection-sha256 "${NT_SHA}" \
  --output "${STAGE}/fit" \
  > "${STAGE}/validation/fit.stdout.json" \
  2> "${STAGE}/validation/fit.stderr.txt"

"${PYTHON}" - "${STAGE}/fit" "${TASKSPEC_SHA}" <<'PY'
import json
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_complete_dna_language_seeded_oof_fit"
    or receipt["candidate_count"] != 11
    or receipt["prediction_rows"] != 113630
    or receipt["task_spec_artifacts_sha256"] != sys.argv[2]
    or receipt["bootstrap_draw_contract"] != "exact_frozen_enformer_sei_offsets"
    or receipt["shared_control_prediction_policy"] != "reuse_bit_exact_frozen_enformer_sei_oof"
    or receipt["shared_control_source_artifacts_sha256"] != "9306d6352163a4ccabf06406709d3545b0a73cca1b47e07ab817bcae874a829e"
    or receipt["shared_control_oof_rows"] != 10330
    or receipt["shared_control_rows_copied_without_numerical_transformation"] is not True
    or receipt["shared_control_refit_performed"] is not False
    or receipt["shared_control_counted_as_new_evidence"] is not False
    or receipt["shared_control_seeds_counted_as_independent_evidence"] is not False
    or receipt["bootstrap_salt_offsets"] != {"inner": 10000, "outer": 20000, "training_mean": 30000}
    or receipt["features_imputed"]
    or receipt["obsolete_6kb_folds_reused"]
    or receipt["metrics_calculated"]
    or receipt["sealed_assets_read"]
):
    raise SystemExit("DNA-language fit firewall differs")
freeze_tree(root, {
    "artifact_class": "gse281364_dna_language_seeded_head_fit",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "models": ["caduceus", "dnabert2", "hyenadna", "nucleotide_transformer"],
    "candidate_count": 11,
    "prediction_rows": 113630,
    "task_spec_artifacts_sha256": sys.argv[2],
    "features_imputed": False,
    "obsolete_6kb_folds_reused": False,
    "metrics_calculated": False,
    "sealed_assets_read": False,
    "status": "pass_complete_dna_language_seeded_oof_fit",
})
verify_frozen_tree(root)
PY
FIT_SHA=$(sha256sum "${STAGE}/fit/ARTIFACTS.json" | awk '{print $1}')

# Refuse family evaluation unless the shared control reproduces its frozen
# Enformer/Sei OOF rows exactly, including prediction text and row ordering.
"${PYTHON}" "${ROOT}/scripts/verify_gse281364_shared_control_identity.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --candidate-fit-tree "${STAGE}/fit" \
  --candidate-fit-sha256 "${FIT_SHA}" \
  --output "${STAGE}/shared_control_identity" \
  > "${STAGE}/validation/shared_control_identity.stdout.json" \
  2> "${STAGE}/validation/shared_control_identity.stderr.txt"

"${PYTHON}" - "${STAGE}/shared_control_identity" "${FIT_SHA}" <<'PY'
import json
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_bit_identical_shared_allele_control_oof_rows"
    or receipt["rows"] != 10330
    or receipt["candidate_fit_artifacts_sha256"] != sys.argv[2]
    or receipt["bit_identical"] is not True
    or receipt["comparison_performed_before_family_evaluation"] is not True
    or receipt["sealed_assets_read"]
):
    raise SystemExit("shared-control identity firewall differs")
freeze_tree(root, {
    "artifact_class": "gse281364_shared_allele_control_identity",
    "dataset_id": "gse281364",
    "rows": 10330,
    "candidate_fit_artifacts_sha256": sys.argv[2],
    "bit_identical": True,
    "comparison_performed_before_family_evaluation": True,
    "sealed_assets_read": False,
    "status": "pass_bit_identical_shared_allele_control_oof_rows",
})
verify_frozen_tree(root)
PY
SHARED_CONTROL_SHA=$(sha256sum "${STAGE}/shared_control_identity/ARTIFACTS.json" | awk '{print $1}')

"${PYTHON}" "${ROOT}/scripts/evaluate_gse281364_dna_language_seeded_heads.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --task-spec-tree "${STAGE}/task_spec" \
  --task-spec-sha256 "${TASKSPEC_SHA}" \
  --fit-tree "${STAGE}/fit" \
  --fit-sha256 "${FIT_SHA}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json" \
  2> "${STAGE}/validation/evaluation.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum "${SOURCE_FILES[@]}" > "${STAGE}/source.terminal.sha256"
cmp --silent "${STAGE}/source.start.sha256" "${STAGE}/source.terminal.sha256"
sha256sum --check --strict "${STAGE}/source.terminal.sha256" > "${STAGE}/validation/terminal_source_recheck.stdout.txt"

"${PYTHON}" - "${STAGE}" "${NT_SHA}" "${TASKSPEC_SHA}" "${FIT_SHA}" "${SHARED_CONTROL_SHA}" <<'PY'
import json
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
evaluation = json.loads((stage / "evaluation/receipt.json").read_text(encoding="utf-8"))
if (
    evaluation["status"] != "pass_independent_dna_language_seeded_development_evaluation"
    or evaluation["cross_family_ranking"]
    or evaluation["overall_leader_reported"]
    or evaluation["features_imputed"]
    or evaluation["obsolete_6kb_folds_reused"]
    or evaluation["shortlist_created"]
    or evaluation["champion_claim"]
    or evaluation["sealed_assets_read"]
):
    raise SystemExit("DNA-language terminal firewall differs")
freeze_tree(stage, {
    "artifact_class": "gse281364_dna_language_seeded_head_campaign",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "models": ["caduceus", "dnabert2", "hyenadna", "nucleotide_transformer"],
    "candidate_count": 11,
    "prediction_rows": 113630,
    "elements": 1033,
    "long_range_blocks": 239,
    "nucleotide_transformer_projection_artifacts_sha256": sys.argv[2],
    "task_spec_artifacts_sha256": sys.argv[3],
    "fit_artifacts_sha256": sys.argv[4],
    "shared_control_identity_artifacts_sha256": sys.argv[5],
    "shared_control_oof_rows_bit_identical_to_enformer_sei": True,
    "features_imputed": False,
    "obsolete_6kb_folds_reused": False,
    "models_reported_separately": True,
    "cross_family_ranking": False,
    "mandatory_task_native_baselines_complete": False,
    "shortlist_created": False,
    "champion_claim": False,
    "sealed_evaluation": False,
    "sealed_assets_read": False,
    "global_frozen_census_changed": False,
    "status": "passed_dna_language_seeded_development_campaign",
})
verify_frozen_tree(stage)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Completed DNA-language seeded development campaign: %s\n' "${TARGET}"
