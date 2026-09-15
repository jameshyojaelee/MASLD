#!/usr/bin/env bash
#SBATCH --job-name=model-probe-014
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --gres=gpu:b6k:1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# First GPU validation of the scvi-tools 1.4.2 inductive runtime, on B6K.
#
# This is NOT a re-validation across a change.  The runtime's only existing
# validation is CPU-only: executions/environments/scvi-tools-1.4.2-inductive-
# portable-21079606/validation/portable-inductive-probe/probe.json records
# "cuda_available": false.  Its torch is 2.13.0+cu130, whose arch list is
# sm_75/80/86/90/100/120 with no sm_89, so it cannot initialise CUDA on an L40S
# and has therefore never executed on any GPU.  B6K is its first.
#
# The probe script and the inductive unit tests are run UNMODIFIED and bound by
# SHA-256, so the only difference from the recorded CPU run is the hardware.
# scvi-tools defaults to accelerator="auto", so a visible CUDA device is used.
#
# The state digests will NOT match the CPU run and are not expected to: different
# kernels give different weights.  What must hold on GPU is the same set of
# invariants the CPU probe asserted -- scvi_state_unchanged, scanvi_state_unchanged,
# direct inference for both models, no query training, no query label exposure,
# no query adaptation, normalised probabilities.
#
# Submitted directly rather than through the GPU dependency dispatcher, which
# hardcodes gpu:l40s:1.  Job ids from this wrapper never enter the dispatcher
# roster, so job 21100574 and the four inherited lanes are unaffected.
#
# Inherited job 21100293 (model-work-151) is deliberately LEFT PENDING and is not
# cancelled, requeued or modified: it is a lane root with nine dependents chained
# behind it by afterany, and letting it terminate on its own inline
# torch.cuda.is_available() guard satisfies those dependencies unambiguously,
# whereas cancelling a pending root risks stranding the chain.
#
# Synthetic 40-gene / 120-row fixture only.  Reads no project data, no development
# outcome, no held-back outcome, and calculates no benchmark metric.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-portable-21079606"
RUNTIME_ARTIFACTS_SHA256=d885ada24fdad9f6e3d57a92ff870ca8ff05f07701ac6b0762f55a8cb19db61a
PROBE="${ROOT}/scripts/probe_scvi_scanvi_1_4_2_inductive.py"
PROBE_SHA256=dccd8ae73483ae3c7e0047a01fd1d1301b42ccba9a11ebac149814b3d7df4362
CPU_PROBE="${RUNTIME}/validation/portable-inductive-probe/probe.json"

TARGET="${ROOT}/executions/scvi-scanvi-1-4-2-b6k-first-gpu-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite validation artifact: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}"

cd "${ROOT}"
sha256sum --check --strict <<EOF
${RUNTIME_ARTIFACTS_SHA256}  ${RUNTIME}/ARTIFACTS.json
${PROBE_SHA256}  ${PROBE}
EOF

nvidia-smi > "${STAGE}/nvidia-smi.txt" 2>&1 || true

module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH=${RUNTIME}/site-packages:${ROOT}:${ROOT}/src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Layer guard: prove the intended torch layer loaded and that a B6K device is live.
python - "${STAGE}" <<'PYA'
import json, os, sys
from pathlib import Path
import torch

stage = Path(sys.argv[1])
expected = "scvi-tools-1.4.2-inductive-portable-21079606"
resolved = os.path.realpath(torch.__file__)
guard = {
    "arch_list": list(torch.cuda.get_arch_list()),
    "cuda_available": torch.cuda.is_available(),
    "cuda_runtime": torch.version.cuda,
    "expected_layer": expected,
    "layer_guard_passed": expected in resolved,
    "torch_file": resolved,
    "torch_version": torch.__version__,
}
if torch.cuda.is_available():
    guard["device_name"] = torch.cuda.get_device_name(0)
    guard["device_capability"] = list(torch.cuda.get_device_capability(0))
    guard["device_count"] = torch.cuda.device_count()
(stage / "layer_guard.json").write_text(json.dumps(guard, indent=2, sort_keys=True) + "\n")
print(json.dumps(guard, indent=2, sort_keys=True))
if not guard["layer_guard_passed"]:
    raise SystemExit(f"layer guard failed: {resolved}")
if torch.__version__ != "2.13.0+cu130":
    raise SystemExit(f"torch is {torch.__version__}")
if not guard["cuda_available"] or guard.get("device_count") != 1:
    raise SystemExit("exactly one visible CUDA device is required")
if guard["device_capability"][0] != 12:
    raise SystemExit(f"expected a Blackwell sm_12x device, got {guard['device_name']}")
if "sm_89" in guard["arch_list"]:
    raise SystemExit("unexpected sm_89 SASS; this build was recorded as L40S-incapable")
PYA

python -m unittest -v tests.unit.test_scvi_scanvi_inductive \
  > "${STAGE}/unit-tests.txt" 2>&1

python "${PROBE}" --output "${STAGE}/b6k-inductive-probe" \
  > "${STAGE}/b6k-inductive-probe.stdout.json" \
  2> "${STAGE}/b6k-inductive-probe.stderr.txt"

# Cross-PROCESS determinism: the probe fixes scvi.settings.seed = 20260824 and a
# seeded numpy generator, so a second independent process must reproduce the same
# module-state digests.  Run in a separate process, not a loop, so nothing carries
# over in memory.
python "${PROBE}" --output "${STAGE}/b6k-inductive-probe-repeat" \
  > "${STAGE}/b6k-inductive-probe-repeat.stdout.json" \
  2> "${STAGE}/b6k-inductive-probe-repeat.stderr.txt"

# Compare invariants (not weights) against the recorded CPU-only probe.
python - "${STAGE}" "${CPU_PROBE}" <<'PYB'
import json, sys
from pathlib import Path

stage, cpu_probe = Path(sys.argv[1]), Path(sys.argv[2])
gpu = json.loads((stage / "b6k-inductive-probe" / "probe.json").read_text())
rep = json.loads((stage / "b6k-inductive-probe-repeat" / "probe.json").read_text())
cpu = json.loads(cpu_probe.read_text())

INVARIANTS = (
    "direct_scanvi_inference",
    "direct_scvi_inference",
    "knn_probabilities_normalized",
    "scanvi_probabilities_normalized",
    "scanvi_state_unchanged",
    "scvi_state_unchanged",
)
MUST_BE_FALSE = (
    "project_data_read",
    "query_adaptation_executed",
    "query_labels_exposed_to_model",
    "query_training_executed",
    "sealed_outcomes_read",
)
failures = []
if gpu.get("status") != "pass":
    failures.append(f"probe status {gpu.get('status')}")
if gpu.get("cuda_available") is not True:
    failures.append("probe did not run on a CUDA device")
if cpu.get("cuda_available") is not False:
    failures.append("reference probe was not the CPU-only record")
for key in INVARIANTS:
    if gpu.get(key) is not True:
        failures.append(f"{key} is not True on GPU")
    if cpu.get(key) is not True:
        failures.append(f"{key} is not True in the CPU record")
for key in MUST_BE_FALSE:
    if gpu.get(key) is not False:
        failures.append(f"{key} is not False on GPU")
for key in ("scvi_tools", "torch", "genes", "train_rows", "query_rows"):
    if gpu.get(key) != cpu.get(key):
        failures.append(f"{key} differs: gpu={gpu.get(key)} cpu={cpu.get(key)}")

# Cross-process determinism, run on the same device in a fresh interpreter.
determinism = {
    "scvi_state_sha256_matches": gpu.get("scvi_state_sha256") == rep.get("scvi_state_sha256"),
    "scanvi_state_sha256_matches": gpu.get("scanvi_state_sha256") == rep.get("scanvi_state_sha256"),
    "repeat_status": rep.get("status"),
    "repeat_cuda_available": rep.get("cuda_available"),
}
cross_process_ok = all(
    determinism[k] for k in ("scvi_state_sha256_matches", "scanvi_state_sha256_matches")
) and rep.get("status") == "pass" and rep.get("cuda_available") is True
if not cross_process_ok:
    failures.append("cross-process determinism failed")

summary = {
    "schema_version": "masld-bench-scvi-scanvi-b6k-first-gpu-validation-v1",
    "cross_process_determinism": determinism,
    "cross_process_determinism_passed": cross_process_ok,
    "cross_process_attribution_note": (
        "The probe does NOT enable torch.use_deterministic_algorithms, so a "
        "divergence here would be a property of the scvi-tools 1.4.2 runtime on "
        "GPU, not a B6K defect. Attribute any failure to the runtime and report "
        "the differing digests rather than the hardware."
    ),
    "repeat_probe_scvi_state_sha256": rep.get("scvi_state_sha256"),
    "repeat_probe_scanvi_state_sha256": rep.get("scanvi_state_sha256"),
    "claim": "FIRST GPU EXECUTION and FIRST GPU VALIDATION of this runtime. Behaviour is being established from nothing, not re-established across a change; there is no prior GPU result to diff against.",
    "cpu_reference_probe": cpu_probe.as_posix(),
    "cpu_reference_cuda_available": cpu.get("cuda_available"),
    "gpu_probe_cuda_available": gpu.get("cuda_available"),
    "invariants_checked": list(INVARIANTS),
    "must_be_false_checked": list(MUST_BE_FALSE),
    "state_digests_compared_across_devices": False,
    "state_digest_note": (
        "scvi_state_sha256 and scanvi_state_sha256 are NOT compared with the CPU "
        "record: different kernels give different weights. Only within-run "
        "invariance (query did not mutate the frozen state) is asserted."
    ),
    "gpu_scvi_state_sha256": gpu.get("scvi_state_sha256"),
    "gpu_scanvi_state_sha256": gpu.get("scanvi_state_sha256"),
    "cpu_scvi_state_sha256": cpu.get("scvi_state_sha256"),
    "cpu_scanvi_state_sha256": cpu.get("scanvi_state_sha256"),
    "l40s_incapability_evidence": (
        "torch 2.13.0+cu130 arch list has no sm_89; b6k admission cell 21106572 "
        "recorded this build reporting cuda_available=false on driver 550.107.02. "
        "CAVEAT: that cell used the sibling scvi-tools-1.5.0.post1-portable-21064579 "
        "build, confirmed byte-identical by torch __version__ and git_version "
        "cf30153c4c131c8164ee7798e5022d810682e2cb. The 1.4.2-inductive tree was NOT "
        "itself executed on an L40S. This is an inference, not a verified fact."
    ),
    "preexisting_queue_defect": (
        "PRE-EXISTING, not introduced by this work: bundle "
        "scvi-scanvi-50k-study-held-3seed-v1 was queued as job 21100293 on "
        "--gres=gpu:l40s:1 with a runtime that has neither sm_89 SASS nor any GPU "
        "validation of any kind (its only validation records "
        "cuda_visible_on_cpu_node=false). Queue membership was therefore never "
        "proof of executability."
    ),
    "authoritative_result_set": (
        "DECIDED BEFORE EITHER RESULT EXISTS: if job 21100293 ever produces a "
        "result set for these 15 logical tasks, the B6K bundle's output is "
        "authoritative, because it is the only one produced on hardware where the "
        "runtime has a GPU validation. 21100293's output would be admissible only "
        "after an L40S GPU validation of this runtime, which does not exist."
    ),
    "inherited_job_disposition": (
        "Job 21100293 (model-work-151, bundle scvi-scanvi-50k-study-held-3seed-v1) "
        "is left PENDING and untouched. It is a GPU lane root with nine afterany "
        "dependents; letting it terminate on its own torch.cuda.is_available() "
        "guard satisfies those dependencies, whereas cancelling a pending root "
        "risks stranding the chain."
    ),
    "status": "fail" if failures else "pass",
}
if failures:
    summary["failures"] = failures
(stage / "validation.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary, indent=2, sort_keys=True))
if failures:
    raise SystemExit("B6K first-GPU validation failed: " + "; ".join(failures))
PYB

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" "${PROBE_SHA256}" <<'PYC'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "scvi_scanvi_b6k_first_gpu_validation",
        "model_ids": ["scvi", "scanvi"],
        "probe_source_sha256": sys.argv[3],
        "runtime_artifacts_sha256": sys.argv[2],
        "validation_sha256": sha256_file(root / "validation.json"),
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PYC

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'scvi/scanvi 1.4.2 B6K first-GPU validation complete: %s\n' "${TARGET}"
