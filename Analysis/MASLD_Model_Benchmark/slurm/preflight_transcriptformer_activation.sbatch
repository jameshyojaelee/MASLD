#!/usr/bin/env bash
#SBATCH --job-name=tf-preflight
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
runtime="${root}/executions/environments/transcriptformer-torch2.5.1-21070738"
runtime_artifacts_sha256="f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb"
python="${runtime}/env/bin/python"
checkpoint="${root}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens"
source="${root}/executions/resource-atlas-smoke-20260821T221938Z-20984634/subset/resource_atlas_geneformer_smoke_1000.h5ad"
vocab="${checkpoint}/vocabs/homo_sapiens_gene.h5"
target="${root}/executions/transcriptformer-activation-preflight-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer activation preflight.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}/validation" "${stage}/cache"
export XDG_CACHE_HOME="${stage}/cache"
export MPLCONFIGDIR="${stage}/cache/matplotlib"
printf '%s  %s\n' "${runtime_artifacts_sha256}" "${runtime}/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${root}"
PYTHONPATH="${root}:${root}/src" "${python}" -m unittest \
  tests.unit.test_build_transcriptformer_unlabeled_fixture \
  > "${stage}/validation/unittest.txt" 2>&1
PYTHONPATH="${root}:${root}/src" "${python}" \
  scripts/build_transcriptformer_unlabeled_fixture.py \
  --source "${source}" \
  --vocab "${vocab}" \
  --output "${stage}/fixture" \
  > "${stage}/validation/fixture.stdout.json"
PYTHONPATH="${root}:${root}/src" "${python}" \
  scripts/preflight_transcriptformer_cpu.py \
  --checkpoint-root "${checkpoint}" \
  --fixture-root "${stage}/fixture" \
  --output "${stage}/validation/cpu_preflight.json" \
  > "${stage}/validation/cpu_preflight.stdout.json"
"${python}" -m pip freeze --all > "${stage}/validation/pip-freeze.txt"
sha256sum \
  config/adaptation/transcriptformer.toml \
  scripts/build_transcriptformer_unlabeled_fixture.py \
  scripts/preflight_transcriptformer_cpu.py \
  scripts/transcriptformer_extract_unlabeled_embeddings.py \
  tests/unit/test_build_transcriptformer_unlabeled_fixture.py \
  slurm/preflight_transcriptformer_activation.sbatch \
  "${runtime}/ARTIFACTS.json" \
  > "${stage}/source.sha256"
PYTHONPATH="${root}/src" "${python}" - "${stage}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
fixture = json.loads((root / "fixture/fixture_receipt.json").read_text())
preflight = json.loads((root / "validation/cpu_preflight.json").read_text())
if (
    fixture["status"] != "pass_outcome_blind_fixture"
    or fixture["evaluation_label_columns_read"]
    or fixture["sealed_outcomes_read"]
    or preflight["status"] != "pass_cpu_preflight"
    or not preflight["checkpoint_strict_load"]
    or preflight["model_parameter_count"] != 429870637
    or preflight["evaluation_label_columns_read"]
    or preflight["sealed_outcomes_read"]
):
    raise SystemExit("TranscriptFormer CPU preflight firewall differs")
freeze_tree(root, {
    "artifact_class": "transcriptformer_outcome_blind_activation_preflight",
    "checkpoint_strict_load": True,
    "dataset_view_id": "resource_atlas_geneformer_smoke_1000_v1",
    "development_rows": 1000,
    "donors": 102,
    "evaluation_labels_read": False,
    "exposure_status": "unknown",
    "gene_order_policies": ["native", "common"],
    "sealed_champion_eligible": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
mv "${stage}" "${target}"
printf 'TranscriptFormer activation CPU preflight: %s\n' "${target}"
sha256sum "${target}/ARTIFACTS.json"
