#!/usr/bin/env bash
#SBATCH --job-name=model-check-998
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
TARGET="${ROOT}/executions/gse281367-label-free-query-atac-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite query-ATAC validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "01c2c38d7583fa16c06444b8870e6bccd181f4e2d1ad57f05749da820cdbe51b config/campaigns/gse281367_label_free_query_atac_materialization_20260825.json" \
  "164847977304d0d545e6448ed8a7371bce4a353ca0d64ff52273f0557e985776 scripts/materialize_gse281367_label_free_query_atac.py" \
  "88e957fc3ff3ab4dcb70eed25c94d4d8601369e352655587ff53c9da5e4daefa tests/unit/test_materialize_gse281367_label_free_query_atac.py" \
  "c7a5eaf49784aa76b13a1cd9edaa94220775c7378ab9272fb0c2bdcd1d3a9873 executions/gse281367-label-free-atac-exchange-axis-20260825/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_materialize_gse281367_label_free_query_atac \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/campaigns/gse281367_label_free_query_atac_materialization_20260825.json \
  scripts/materialize_gse281367_label_free_query_atac.py \
  tests/unit/test_materialize_gse281367_label_free_query_atac.py \
  slurm/validate_gse281367_label_free_query_atac_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" "${SLURM_JOB_ID}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
stderr = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if "Ran 5 tests" not in stderr or not stderr.rstrip().endswith("OK"):
    raise SystemExit("query-ATAC materialization tests did not pass exactly")
freeze_tree(
    root,
    {
        "artifact_class": "gse281367_label_free_query_atac_validation",
        "slurm_job_id": sys.argv[2],
        "tests": 5,
        "scored_target_contigs_excluded": True,
        "condition_values_read": False,
        "evaluator_outcome_artifact_read": False,
        "metrics_calculated": False,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Query-ATAC validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
