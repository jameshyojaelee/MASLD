#!/usr/bin/env bash
#SBATCH --job-name=data-prep-264
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
PREFLIGHT="${ROOT}/executions/model-check-260-21096725"
TILES="${ROOT}/executions/corgi-tile-smoke-subset-21068454"
TARGET="${ROOT}/executions/data-prep-264-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi FiLM-head fixture v2: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for pair in \
  "01de778c17c67d92de9a2e242350e17443dbfe76b7c5a5f8bc447eae5408cd3c ${RUNTIME}/ARTIFACTS.json" \
  "eace72d1b831a4834d692134749bb82093c29930fb3a0258ae667b8acb840959 ${IMAGE}" \
  "378a21b2ccb7b6e9cae865cd11aef61309f8a3f6e5f49bedbe640b5aca4a82a7 ${PREFLIGHT}/ARTIFACTS.json" \
  "2617f3a85a9ced8b0f7cb9efe3819b1e57177629221a2013e39fc7c96e187f5b ${TILES}/ARTIFACTS.json" \
  "d131cb83b3d0391364515c905e1dd4d1c6f7447c30024aa041c5ae46ec2cae70 ${ROOT}/scripts/prepare_corgi_regular_film_head_fixture.py" \
  "ba6dc44e1e511e5c94217efaf3a3abac67d1baa5bfca3a7df8c9ff45c57c4655 ${ROOT}/tests/unit/test_prepare_corgi_regular_film_head_fixture.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module purge
module load singularity/4.2.2
COMMON=(--cleanenv --containall --net --network none --bind "${ROOT}:${ROOT}:ro" --bind "${STAGE}:${STAGE}:rw" --pwd "${ROOT}" --env "PYTHONPATH=${ROOT}" --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 --env PYTHONHASHSEED=20260825)
singularity exec "${COMMON[@]}" "${IMAGE}" python -m unittest -v \
  tests.unit.test_prepare_corgi_regular_film_head_fixture \
  > "${STAGE}/validation/tests.stdout.txt" 2> "${STAGE}/validation/tests.stderr.txt"
singularity exec "${COMMON[@]}" "${IMAGE}" python \
  scripts/prepare_corgi_regular_film_head_fixture.py \
  --contract-preflight "${PREFLIGHT}" --tile-roster "${TILES}" \
  --mapper "0=${ROOT}/executions/corgi-gse296875-mapper-fold0-21066334" \
  --mapper "1=${ROOT}/executions/corgi-gse296875-mapper-fold1-21066335" \
  --mapper "2=${ROOT}/executions/corgi-gse296875-mapper-fold2-21066337" \
  --mapper "3=${ROOT}/executions/corgi-gse296875-mapper-fold3-21066336" \
  --mapper "4=${ROOT}/executions/corgi-gse296875-mapper-fold4-21066338" \
  --output "${STAGE}/fixture" > "${STAGE}/validation/prepare.stdout.json"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --pwd "${ROOT}" "${IMAGE}" \
  python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
cd "${ROOT}"
sha256sum scripts/prepare_corgi_regular_film_head_fixture.py \
  tests/unit/test_prepare_corgi_regular_film_head_fixture.py \
  slurm/prepare_corgi_regular_film_head_fixture_cpu_v2.sbatch > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "fixture/receipt.json").read_text())
closed = (
    "atac_bigwig_manifest_or_signal_read",
    "development_atac_outcomes_read",
    "histology_or_disease_labels_read",
    "test_or_sealed_features_or_outcomes_read",
    "model_checkpoint_loaded",
    "model_fit_performed",
    "model_prediction_performed",
    "benchmark_metrics_computed",
    "training_execution_authorized",
    "prediction_execution_authorized",
    "evaluation_execution_authorized",
)
if (
    receipt["status"] != "pass_outcome_blind_crossed_fold_fixture"
    or len(receipt["folds"]) != 5
    or receipt["contexts"] != 195
    or receipt["donors"] != 39
    or receipt["context_width"] != 2891
    or any(receipt[field] for field in closed)
    or any(
        fold["fit_held_donor_overlap"] != 0
        or fold["fit_held_genomic_fold_overlap"] != 0
        for fold in receipt["folds"]
    )
):
    raise SystemExit("Corgi FiLM-head fixture receipt differs")
freeze_tree(root, {
    "artifact_class": "corgi_regular_film_plus_head_outcome_blind_fixture_v2",
    "model_id": "corgi_regular", "task_id": "rna_conditioned_atac",
    "dataset_id": "gse296875", "lineage_id": "hepatocyte",
    "adaptation_rung": "film_plus_head", "biological_unit": "donor",
    "genomic_block_unit": "chromosome", "folds": 5,
    "donors": receipt["donors"], "contexts": receipt["contexts"],
    "schedule_rows": receipt["schedule_rows"], "tests": 4,
    "supersedes_failed_attempt_job_id": "21096816",
    "failed_attempt_disposition": "one_character_fold4_mapper_sha_typo_before_mapper_contents",
    "development_atac_outcomes_read": False,
    "histology_or_disease_labels_read": False,
    "test_or_sealed_features_or_outcomes_read": False,
    "model_checkpoint_loaded": False, "model_fit_performed": False,
    "model_prediction_performed": False, "benchmark_metrics_computed": False,
    "training_execution_authorized": False, "prediction_execution_authorized": False,
    "evaluation_execution_authorized": False, "global_census_modified": False,
    "global_promotion_gate_modified": False, "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi FiLM-plus-head outcome-blind fixture v2: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
