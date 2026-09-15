#!/usr/bin/env bash
#SBATCH --job-name=tf-queue-check
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
queue="${root}/config/campaigns/gpu_bundle_queue/042-transcriptformer-activation-retry.json"
target="${root}/executions/transcriptformer-gpu-queue-validation-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer queue validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}/validation"
module load python/3.11.5-GCCcore-13.2.0
cd "${root}"
PYTHONPATH="${root}/src:${root}/scripts" python - \
  "${queue}" "${stage}/validation/receipt.json" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import sys

from gpu_bundle_dispatcher import load_item
from masld_bench.artifacts import verify_frozen_tree

queue_path = Path(sys.argv[1]).resolve(strict=True)
receipt_path = Path(sys.argv[2])
root = Path.cwd().resolve(strict=True)
item = load_item(queue_path, root)
if (
    item["bundle_id"] != "transcriptformer-activation-retry1"
    or item["enabled"]
    or item["priority"] != 42
    or item["logical_tasks"] != 2
    or item["family"] != "cell_foundation"
    or not item["_ready"]
):
    raise SystemExit("TranscriptFormer disabled queue item differs")
wrapper = item["_wrapper"].read_text()
for required in (
    "#SBATCH --job-name=model-training",
    "#SBATCH --partition=gpu",
    "#SBATCH --account=nslab",
    "#SBATCH --qos=nslab",
    "#SBATCH --gres=gpu:l40s:1",
):
    if wrapper.count(required) != 1:
        raise SystemExit(f"TranscriptFormer GPU wrapper differs: {required}")
for required in ("evaluation_label_columns_read", "sealed_outcomes_read"):
    if required not in wrapper:
        raise SystemExit(f"TranscriptFormer GPU firewall differs: {required}")
for forbidden in ("--qos=innovation", "#SBATCH --array", "#SBATCH --job-name=masld"):
    if forbidden in wrapper:
        raise SystemExit(f"forbidden GPU wrapper value: {forbidden}")
expected_required = {
    "executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json":
        "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb",
    "executions/transcriptformer-activation-preflight-21070971/ARTIFACTS.json":
        "2402298e4492a487c063962b4b79d95fe15e871670d972a356a8b24b0a6a588c",
}
if set(item["required_paths"]) != set(expected_required):
    raise SystemExit("TranscriptFormer required artifacts differ")
for relative, expected_sha256 in expected_required.items():
    artifact = (root / relative).resolve(strict=True)
    if sha256(artifact.read_bytes()).hexdigest() != expected_sha256:
        raise SystemExit("TranscriptFormer required artifact hash differs")
    if "activation-preflight" in relative:
        verify_frozen_tree(artifact.parent)
runtime = json.loads(
    (root / "executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json").read_text()
)
if runtime.get("status") != "pass_cpu_runtime" or runtime.get("transcriptformer_version") != "0.6.1":
    raise SystemExit("TranscriptFormer runtime receipt differs")
receipt = {
    "bundle_id": item["bundle_id"],
    "enabled": False,
    "gpu_submitted": False,
    "logical_tasks": item["logical_tasks"],
    "outcomes_read": False,
    "queue_sha256": sha256(queue_path.read_bytes()).hexdigest(),
    "schema_version": "masld-bench-transcriptformer-gpu-queue-validation-v1",
    "status": "pass_disabled_queue_item",
    "wrapper_sha256": item["wrapper_sha256"],
}
receipt_path.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
PY
python -m pip freeze --all > "${stage}/validation/pip-freeze.txt"
sha256sum \
  config/campaigns/gpu_bundle_queue/042-transcriptformer-activation-retry.json \
  scripts/gpu_bundle_dispatcher.py \
  slurm/run_transcriptformer_activation_retry.sbatch \
  slurm/validate_transcriptformer_gpu_queue.sbatch \
  > "${stage}/source.sha256"
PYTHONPATH="${root}/src" python - "${stage}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/receipt.json").read_text())
if receipt["status"] != "pass_disabled_queue_item" or receipt["enabled"] or receipt["gpu_submitted"]:
    raise SystemExit("TranscriptFormer queue validation receipt differs")
freeze_tree(root, {
    "artifact_class": "transcriptformer_gpu_queue_validation",
    "bundle_id": receipt["bundle_id"],
    "enabled": False,
    "gpu_submitted": False,
    "outcomes_read": False,
    "status": "passed",
})
PY
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
