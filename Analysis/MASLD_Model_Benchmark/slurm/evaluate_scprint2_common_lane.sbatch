#!/usr/bin/env bash
#SBATCH --job-name=model-evaluation-302
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
python="${root}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
embedding_artifact="${root}/executions/model-work-302-21075963"
embedding_artifact_sha256="3c57bcf5f9c821e395966809033b4319971432102714c0051f2c8d449e654183"
atlas="${root}/executions/resource-atlas-smoke-20260821T221938Z-20984634/subset/resource_atlas_geneformer_smoke_1000.h5ad"
atlas_sha256="d59c85b4fb0a1352a427de3b2a587c30aac12891a8d753e10d4442fc719c985f"
row_id_source="${root}/executions/scprint2-no-neighbor-fixture-21075645/fixture/aligned_counts.h5ad"
row_id_source_sha256="951b2615c437c86991bf039b885c453c51d9af3b3d39a111cf33a914d8489635"
target="${root}/executions/model-evaluation-302-${SLURM_JOB_ID}"
stage="${target}.staging"
failed="${target}.failed"

preserve_failed_stage() {
  if [[ -d "${stage}" && ! -e "${failed}" ]]; then
    mv "${stage}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${target}" || -e "${stage}" || -e "${failed}" ]]; then
  printf 'Refusing to overwrite model-evaluation-302 artifact.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}/evaluation" "${stage}/validation"
printf '%s  %s\n' "${embedding_artifact_sha256}" \
  "${embedding_artifact}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${atlas_sha256}" "${atlas}" | sha256sum --check --strict
printf '%s  %s\n' "${row_id_source_sha256}" "${row_id_source}" | sha256sum --check --strict

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PYTHONPATH="${root}/src"
cd "${root}"

representations=(
  common
  random_checkpoint_max
  released_selection_default_no_neighbor
)
for representation in "${representations[@]}"; do
  "${python}" scripts/evaluate_transcriptformer_common_lane.py \
    --embeddings "${embedding_artifact}/embeddings/${representation}_embeddings.npz" \
    --atlas "${atlas}" \
    --baseline-root "${root}/executions/folds" \
    --bootstrap-resamples 10000 \
    --model-id scprint2_small_v2 \
    --representation-id "${representation}" \
    --expected-width 424 \
    --row-id-source "${row_id_source}" \
    --schema-version masld-bench-scprint2-common-lane-evaluation-v1 \
    --activation-fold-namespace scprint2-resource-atlas-smoke-no-neighbor-v1 \
    --output "${stage}/evaluation/${representation}" \
    > "${stage}/validation/${representation}.stdout.json" \
    2> "${stage}/validation/${representation}.stderr.txt"
done

"${python}" -m pip freeze --all > "${stage}/validation/pip-freeze.txt"
sha256sum \
  scripts/evaluate_transcriptformer_common_lane.py \
  slurm/evaluate_scprint2_common_lane.sbatch \
  "${embedding_artifact}/ARTIFACTS.json" \
  "${atlas}" \
  "${row_id_source}" \
  > "${stage}/source.sha256"

"${python}" - "${stage}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
representations = {
    "common", "random_checkpoint_max", "released_selection_default_no_neighbor"
}
for representation in representations:
    evaluation = json.loads(
        (root / "evaluation" / representation / "evaluation.json").read_text()
    )
    if (
        evaluation["status"] != "pass_development_smoke"
        or evaluation["model_id"] != "scprint2_small_v2"
        or evaluation["representation_id"] != representation
        or evaluation["development_rows"] != 1000
        or evaluation["development_donors"] != 102
        or evaluation["sealed_outcomes_read"]
        or evaluation["sealed_champion_eligible"]
        or set(evaluation["heads"]) != {"linear", "two_layer_mlp"}
    ):
        raise SystemExit("scPRINT-2 common-head evaluation firewall differs")
freeze_tree(root, {
    "artifact_class": "scprint2_no_neighbor_development_common_lane_evaluation",
    "development_only": True,
    "embedding_artifacts_sha256": "3c57bcf5f9c821e395966809033b4319971432102714c0051f2c8d449e654183",
    "evaluation_labels_read": ["broad_label"],
    "exposure_status": "unknown",
    "model_id": "scprint2_small_v2",
    "neighbor_tensors_used": False,
    "representations": sorted(representations),
    "sealed_champion_eligible": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
