#!/usr/bin/env bash
#SBATCH --job-name=model-evaluation
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
python="${root}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
embedding_artifact="${root}/executions/model-training-21071475"
embedding_artifact_sha256="375f38c25140746a77d13dffc7a90fed47e6597dc0537578bb311d9807043007"
atlas="${root}/executions/resource-atlas-smoke-20260821T221938Z-20984634/subset/resource_atlas_geneformer_smoke_1000.h5ad"
atlas_sha256="d59c85b4fb0a1352a427de3b2a587c30aac12891a8d753e10d4442fc719c985f"
target="${root}/executions/model-evaluation-${SLURM_JOB_ID}"
stage="${target}.staging"
failed="${target}.failed"

preserve_failed_stage() {
  if [[ -d "${stage}" && ! -e "${failed}" ]]; then
    mv "${stage}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${target}" || -e "${stage}" || -e "${failed}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer evaluation.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}/validation"
printf '%s  %s\n' "${embedding_artifact_sha256}" "${embedding_artifact}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${atlas_sha256}" "${atlas}" | sha256sum --check --strict
cd "${root}"
PYTHONPATH="${root}/src" "${python}" \
  scripts/evaluate_transcriptformer_common_lane.py \
  --embeddings "${embedding_artifact}/embeddings/common_embeddings.npz" \
  --atlas "${atlas}" \
  --baseline-root "${root}/executions/folds" \
  --bootstrap-resamples 10000 \
  --output "${stage}/evaluation" \
  > "${stage}/validation/evaluation.stdout.json" \
  2> "${stage}/validation/evaluation.stderr.txt"
"${python}" -m pip freeze --all > "${stage}/validation/pip-freeze.txt"
sha256sum \
  config/adaptation/common_and_native.toml \
  config/adaptation/transcriptformer.toml \
  config/artifacts/evaluator_rosters/cell_state_mapping.json \
  scripts/evaluate_transcriptformer_common_lane.py \
  slurm/evaluate_transcriptformer_common_lane.sbatch \
  "${embedding_artifact}/ARTIFACTS.json" \
  "${atlas}" \
  > "${stage}/source.sha256"
PYTHONPATH="${root}/src" "${python}" - "${stage}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
evaluation = json.loads((root / "evaluation/evaluation.json").read_text())
if (
    evaluation["status"] != "pass_development_smoke"
    or evaluation["development_rows"] != 1000
    or evaluation["development_donors"] != 102
    or evaluation["sealed_outcomes_read"]
    or evaluation["sealed_champion_eligible"]
    or set(evaluation["heads"]) != {"linear", "two_layer_mlp"}
    or set(evaluation["registered_classical_baselines"]) != {
        "hvg_pca_nearest_centroid", "hvg_pca_knn", "hvg_pca_logistic",
        "hvg_pca_elastic_net", "hvg_pca_linear_svm",
    }
    or len(evaluation["existing_cell_foundation_common_lane"]) != 12
):
    raise SystemExit("TranscriptFormer development evaluation firewall differs")
freeze_tree(root, {
    "artifact_class": "transcriptformer_development_common_lane_evaluation",
    "development_only": True,
    "evaluation_labels_read": ["broad_label"],
    "exposure_status": "unknown",
    "model_id": "transcriptformer_tf_sapiens",
    "sealed_champion_eligible": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
