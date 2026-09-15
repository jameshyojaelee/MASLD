#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-359
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-359-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-359-%j.err

# Freeze the validation file that the v12 CPU control execute wrapper binds by
# SHA-256. Runs the focused unit tests and a --preflight (no run executes), then
# records the binding. Mirrors slurm/validate_cobolt_v9_cpu_baseline_bundle.sbatch.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v12--3972f43112b8ecc4
CANDIDATE_SHA=9b878eb08b3fc8a6012f473694b778a771da1bdd976f59e2966776f0bda3519f
STABILITY=${ROOT}/executions/cobolt-v12-revision-scope-validation-21111670
STABILITY_SHA=0952f39b35c6d9ea437a32742bfacaf10d3a513e4507d6812b5ca6174e132e0a
EXECUTOR=${ROOT}/scripts/run_cobolt_v12_cpu_baseline_bundle.py
EXECUTOR_SHA=bd710456aff1e307854b68f0edd214e02f181656c151d35851c1d291e1f8e1fa
TARGET=${ROOT}/executions/cobolt-v12-cpu-control-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite validation output: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/model-cpu-check-359-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EXECUTOR_SHA}" "${EXECUTOR}" | sha256sum --check --strict
export PYTHONPATH=${ROOT}:${ROOT}/src
python -m py_compile "${EXECUTOR}"
bash -n slurm/run_cobolt_v12_cpu_controls.sbatch
python -m unittest -v tests.unit.test_cobolt_v12_cpu_baseline_bundle \
  >"${STAGE}/tests.stdout.txt" 2>"${STAGE}/tests.stderr.txt"
python "${EXECUTOR}" \
  --candidate "${CANDIDATE}" \
  --candidate-sha256 "${CANDIDATE_SHA}" \
  --stability "${STABILITY}" \
  --stability-sha256 "${STABILITY_SHA}" \
  --executor-sha256 "${EXECUTOR_SHA}" \
  --workers 12 \
  --allocated-cpus 16 \
  --allocated-memory-gb 200 \
  --preflight >"${STAGE}/preflight.json"
sha256sum \
  "${EXECUTOR}" \
  tests/unit/test_cobolt_v12_cpu_baseline_bundle.py \
  slurm/run_cobolt_v12_cpu_controls.sbatch \
  slurm/validate_cobolt_v12_cpu_controls.sbatch \
  >"${STAGE}/source.sha256"

python - "${STAGE}" "${CANDIDATE_SHA}" "${STABILITY_SHA}" "${EXECUTOR_SHA}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, write_json_exclusive

stage = Path(sys.argv[1]).resolve(strict=True)
pre = json.loads((stage / "preflight.json").read_text(encoding="utf-8"))
if (
    pre.get("status") != "pass"
    or pre.get("runs") != 150
    or pre.get("execution_started") is not False
    or pre.get("sealed_features_opened") is not False
    or pre.get("sealed_labels_opened") is not False
    or pre.get("outcome_records_parsed_or_joined") is not False
    or pre.get("campaign_id") != "v1-rna-atac-cobolt-5seed-smoke-v12"
):
    raise SystemExit(f"CPU control preflight differs: {json.dumps(pre)[:400]}")
write_json_exclusive(stage / "validation_receipt.json", {
    "schema_version": "masld-bench-cobolt-v12-cpu-control-validation-v1",
    "candidate_manifest_sha256": sys.argv[2],
    "stability_manifest_sha256": sys.argv[3],
    "executor_sha256": sys.argv[4],
    "preflight_runs": 150,
    "execution_started": False,
    "outcomes_opened": False,
    "supersedes_for_v11_comparison_only": (
        "executions/cobolt-v9-current-cpu-baseline-bundles/cobolt-v9-cpu-baselines--21088479"
        " is superseded FOR v12 COMPARISON PURPOSES ONLY. The v9 matrix remains intact"
        " as its own frozen record and is neither deleted nor modified."
    ),
    "why_a_re_run_is_required": (
        "run_id = canonical_sha256(identity_payload), which includes campaign_id and"
        " code_lock_sha256. The B6K runtime revision re-keys every run in the campaign,"
        " including these 150 whose own adapter rna_atac_classical.py did not change."
    ),
    "status": "pass",
})
print(freeze_tree(stage, {
    "artifact_class": "cobolt_v12_cpu_control_validation",
    "candidate_manifest_sha256": sys.argv[2],
    "stability_manifest_sha256": sys.argv[3],
    "executor_sha256": sys.argv[4],
}))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Cobolt v12 CPU control validation frozen: %s\n' "${TARGET}"
