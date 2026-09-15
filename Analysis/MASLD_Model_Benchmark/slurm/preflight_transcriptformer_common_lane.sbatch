#!/usr/bin/env bash
#SBATCH --job-name=model-preflight
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
python="${root}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
target="${root}/executions/model-preflight-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer common-head preflight.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}/validation"
cd "${root}"
"${python}" -m py_compile scripts/evaluate_transcriptformer_common_lane.py
PYTHONPATH="${root}/src" "${python}" - <<'PY' > "${stage}/validation/preflight.json"
from pathlib import Path
import json
import numpy as np
import anndata
from scripts.evaluate_transcriptformer_common_lane import (
    BASELINES, EXISTING_CELL_CANDIDATES, HEADS, ROSTER, _fold_index,
    _donor_class_weights, _inner_assignment,
)

root = Path.cwd()
embedding = root / "executions/model-training-21071475/embeddings/common_embeddings.npz"
atlas_path = root / "executions/resource-atlas-smoke-20260821T221938Z-20984634/subset/resource_atlas_geneformer_smoke_1000.h5ad"
bundle = np.load(embedding, allow_pickle=False)
rows = bundle["row_ids"].astype(str)
folds = bundle["outer_folds"].astype(int)
atlas = anndata.read_h5ad(atlas_path, backed="r")
index = {str(row): i for i, row in enumerate(atlas.obs_names)}
donors = [str(atlas.obs.iloc[index[row]]["donor_id"]) for row in rows]
labels = [str(atlas.obs.iloc[index[row]]["broad_label"]) for row in rows]
atlas.file.close()
if bundle["embeddings"].shape != (1000, 2048):
    raise SystemExit("embedding shape differs")
evaluation_folds = np.asarray([_fold_index(donor) for donor in donors])
for donor in set(donors):
    mask = np.asarray([value == donor for value in donors])
    if len(set(folds[mask])) != 1 or len(set(evaluation_folds[mask])) != 1:
        raise SystemExit("a donor spans activation or evaluation folds")
donors_array = np.asarray(donors)
labels_array = np.asarray(labels)
class_index = {label: index for index, label in enumerate(ROSTER)}
targets = np.asarray([class_index[label] for label in labels_array])
for outer_fold in range(5):
    outer_training = np.flatnonzero(evaluation_folds != outer_fold)
    assignment = _inner_assignment(donors_array[outer_training], outer_fold)
    for inner_fold in range(5):
        fit = [
            index for index in outer_training
            if assignment[donors_array[index]] != inner_fold
        ]
        if set(labels_array[fit]) != set(ROSTER):
            raise SystemExit("an inner fit partition lacks a frozen class")
        _donor_class_weights(
            donors_array[fit], targets[fit], require_full_roster=True
        )
for model in BASELINES:
    for suffix in ("predictions.tsv", "class_probabilities.tsv"):
        paths = sorted((root / "executions/folds" / model).glob(f"fold*/adapter_actions/003-predict/{suffix}"))
        if len(paths) != 5:
            raise SystemExit(f"baseline files differ: {model} {suffix}")
for model in EXISTING_CELL_CANDIDATES:
    for head in HEADS:
        for suffix in ("predictions.tsv", "class_probabilities.tsv"):
            paths = sorted((root / "executions/folds").glob(f"{model}*/{head}/fold*/predict/{suffix}"))
            if len(paths) != 5:
                raise SystemExit(f"candidate files differ: {model} {head} {suffix}")
print(json.dumps({
    "status": "pass",
    "rows": len(rows),
    "donors": len(set(donors)),
    "fold_annotation_mismatch_rows": int(np.sum(folds != evaluation_folds)),
    "activation_fold_annotations_used_for_head_fitting": False,
    "registered_baselines": len(BASELINES),
    "existing_candidate_comparators": len(EXISTING_CELL_CANDIDATES) * len(HEADS),
    "sealed_outcomes_read": False,
}, sort_keys=True))
PY
"${python}" -m pip freeze --all > "${stage}/validation/pip-freeze.txt"
sha256sum \
  scripts/evaluate_transcriptformer_common_lane.py \
  slurm/evaluate_transcriptformer_common_lane.sbatch \
  slurm/preflight_transcriptformer_common_lane.sbatch \
  > "${stage}/source.sha256"
chmod -R a-w "${stage}"
mv "${stage}" "${target}"
sha256sum "${target}/validation/preflight.json"
