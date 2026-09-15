#!/usr/bin/env bash
#SBATCH --job-name=model-check-361
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/b6k-admission-record-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"; FAILED="${TARGET}.failed"
preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve EXIT
[[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]] && { echo "refusing overwrite" >&2; exit 1; }
install -d -m 0750 "${STAGE}" "${STAGE}/records" "${STAGE}/validation"

AUTHORIZATION="${ROOT}/config/b6k_hardware_authorization_20260825.json"
ROUTING="${ROOT}/config/b6k_routing_guidance_20260825.json"
FINETUNE="${ROOT}/config/runtimes/gpu_finetune.toml"
INFERENCE="${ROOT}/config/runtimes/gpu_inference.toml"
TOLERANCE="${ROOT}/config/b6k_numeric_tolerance_prespecification_20260825.json"
RECORD_TEST="${ROOT}/tests/unit/test_b6k_admission_records.py"

# Stage identity check. Every input this output file asserts over is pinned here.
{
  printf '%s  %s\n' 244193dde1cfca83b2d2e702be831447bf274de34b0893d5549af15c9eb0ac33 "${AUTHORIZATION}"
  printf '%s  %s\n' b1fe5f652a00cf5c36a2527377fc339646d2f7743d6a2e5e05491f33cc5c4c28 "${ROUTING}"
  printf '%s  %s\n' c3d7820994b8b3f58934701cacb646405e96a3152c8a4b6a85e3d2b321fb73ce "${FINETUNE}"
  printf '%s  %s\n' 686a5d14abfc146b0ee805c7a6e4bdb0201dd07783b5734a448277e49fec1f4c "${INFERENCE}"
  printf '%s  %s\n' a26a6f3fe0d77ef12706038eaede32f222421154c0972e228b5238febeaa722f "${TOLERANCE}"
  printf '%s  %s\n' 7aa1a07c94afc928d6a2211e88dbd73bbe2440844833e14f3beaf81b8a048580 "${RECORD_TEST}"
} | tee "${STAGE}/validation/staged_inputs.sha256" | sha256sum --check --strict

install -m 0640 "${AUTHORIZATION}" "${STAGE}/records/b6k_hardware_authorization_20260825.json"
install -m 0640 "${ROUTING}" "${STAGE}/records/b6k_routing_guidance_20260825.json"
install -m 0640 "${FINETUNE}" "${STAGE}/records/gpu_finetune.toml"
install -m 0640 "${INFERENCE}" "${STAGE}/records/gpu_inference.toml"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"

PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_b6k_admission_records \
  > "${STAGE}/validation/unit_tests.stdout.txt" \
  2> "${STAGE}/validation/unit_tests.stderr.txt"

# The registry requirements must still hold for every runtime file, not just the two
# that changed. This is the test that would catch a broken missing piece invariant.
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.contract.test_registry_files.RegistryFileContractTest.test_runtime_profiles_are_bounded_and_unresolved_runs_block \
  > "${STAGE}/validation/registry_contract.stdout.txt" \
  2> "${STAGE}/validation/registry_contract.stderr.txt"

# Re-hash the dispatcher's five source bindings from the policy and confirm no
# drift. This job must not perturb running job 21100574.
PYTHONPATH="${ROOT}/src:${ROOT}" python - "${ROOT}" "${STAGE}" <<'PY'
import hashlib, json, sys
from pathlib import Path

root, stage = Path(sys.argv[1]), Path(sys.argv[2])
policy = json.loads(
    (root / "config" / "campaigns" / "gpu_dependency_dispatcher_policy_20260825.json").read_text()
)
observed = {}
for label, record in policy["source_bindings"].items():
    path = root / record["path"]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    observed[label] = {"path": record["path"], "expected": record["sha256"], "observed": digest}
    if digest != record["sha256"]:
        raise SystemExit(f"dispatcher source binding drifted: {label}")
payload = {
    "schema_version": "masld-bench-b6k-record-dispatcher-nondrift-v1",
    "dispatcher_job_id": "21100574",
    "dispatcher_files_touched": False,
    "source_bindings": observed,
    "routing_gres_unchanged": policy["routing"]["gres"],
    "campaign_id_unchanged": policy["campaign_id"],
}
(stage / "validation" / "dispatcher_nondrift.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n"
)
print(json.dumps({k: v["observed"] == v["expected"] for k, v in observed.items()}, sort_keys=True))
PY

python3 -m pip freeze --all > "${STAGE}/pip_freeze.txt" 2>/dev/null || true
sha256sum \
  "${ROOT}/slurm/freeze_b6k_admission_records_cpu.sbatch" \
  "${AUTHORIZATION}" "${ROUTING}" "${FINETUNE}" "${INFERENCE}" \
  "${TOLERANCE}" "${RECORD_TEST}" > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" <<'PY'
import sys
from pathlib import Path
from masld_bench.artifacts import freeze_tree

freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "b6k_admission_authorization_and_routing_record",
    "record_only": True,
    "consumed_by_runtime": False,
    "runtime_config_update_route": "in_place_edit_unpinned",
    "dispatcher_files_touched": False,
    "jobs_cancelled_or_requeued": False,
    "gpu_submitted": False,
    "synthetic_inputs_only": True,
    "biological_data_read": False,
    "development_outcomes_read": False,
    "sealed_outcomes_read": False,
    "benchmark_metrics_calculated": False,
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
