#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-admit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-seed20260824-21064489"
SOURCE_ARTIFACTS_SHA256=a0ed38216f0f0acc139cede8c43060af5884db473c888d6e15092f8b71b3e547
AUDIT_SHA256=d5d94436d151a223c08c57ad52dec04d9b8bba37a8bec503780f97bd6b13ffc9
TARGET="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-seed20260824-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite ChromBPNet admission artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
if [[ "${SOURCE_ARTIFACTS_SHA256}" == TO_BE_FILLED* ]]; then
  printf 'Source artifact binding is unresolved.\n' >&2
  exit 1
fi

printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" \
  "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${AUDIT_SHA256}" \
  "${ROOT}/scripts/audit_chrombpnet_prediction_artifact.py" \
  | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
python --version > "${ROOT}/executions/masld-chrombpnet-admit-${SLURM_JOB_ID}.python.txt"
micromamba list --explicit \
  > "${ROOT}/executions/masld-chrombpnet-admit-${SLURM_JOB_ID}.environment.txt"
PYTHONPATH="${ROOT}/src" python \
  "${ROOT}/scripts/audit_chrombpnet_prediction_artifact.py" \
    --artifact "${SOURCE}" \
    --expected-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
    --output "${STAGE}" \
    > "${ROOT}/executions/masld-chrombpnet-admit-${SLURM_JOB_ID}.stdout.txt" \
    2> "${ROOT}/executions/masld-chrombpnet-admit-${SLURM_JOB_ID}.stderr.txt"

mv "${ROOT}/executions/masld-chrombpnet-admit-${SLURM_JOB_ID}.python.txt" \
  "${STAGE}/python.version.txt"
mv "${ROOT}/executions/masld-chrombpnet-admit-${SLURM_JOB_ID}.environment.txt" \
  "${STAGE}/environment.explicit.txt"
mv "${ROOT}/executions/masld-chrombpnet-admit-${SLURM_JOB_ID}.stdout.txt" \
  "${STAGE}/audit.stdout.txt"
mv "${ROOT}/executions/masld-chrombpnet-admit-${SLURM_JOB_ID}.stderr.txt" \
  "${STAGE}/audit.stderr.txt"
sha256sum \
  "${ROOT}/scripts/audit_chrombpnet_prediction_artifact.py" \
  "${ROOT}/slurm/audit_chrombpnet_hepatocyte_donor0_genomic0_seed20260824.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${SOURCE_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

print(
    freeze_tree(
        Path(sys.argv[1]),
        {
            "artifact_class": "chrombpnet_prediction_admission",
            "model_id": "chrombpnet",
            "dataset_id": "gse296875",
            "lineage_id": "hepatocyte",
            "split_id": "donor0_genomic0",
            "seed": 20260824,
            "adaptation_lane": "native_full_head",
            "source_artifacts_sha256": sys.argv[2],
            "primary_biological_score": "nobias_model",
            "assay_qc_score": "full_model",
            "admission_disposition": (
                "admitted_for_independent_development_evaluation_only"
            ),
            "champion_eligible": False,
            "evaluator_outcomes_exposed": False,
            "benchmark_metrics_calculated": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet prediction admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
