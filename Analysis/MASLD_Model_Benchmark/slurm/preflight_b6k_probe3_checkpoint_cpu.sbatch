#!/usr/bin/env bash
#SBATCH --job-name=b6kckpt
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Probe 3 is entirely CPU-resident by construction: it loads with
# map_location="cpu" and checks strict restore semantics, none of which needs a
# device.  Running it here first means a shape or key-structure bug in the
# probe surfaces on a free CPU slot instead of burning a scarce GPU allocation.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH_OVERLAY="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
COBOLT_CKPT="${ROOT}/executions/cobolt-runtime-probe-21060919/state_dict.pt"
COBOLT_SHA256=fa329ed1f4649c1cc375ca042211e09ee74402d4b1ba5df9aaca0a56edbadaeb

TARGET="${ROOT}/executions/b6k-preflight-probe3-checkpoint-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite preflight artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

sha256sum --check --strict <<EOF
${COBOLT_SHA256}  ${COBOLT_CKPT}
30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af  ${BASE}/ARTIFACTS.json
d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3  ${TORCH_OVERLAY}/ARTIFACTS.json
EOF

module load python/3.12.3-GCCcore-13.3.0

for cell in "cu124-overlay-stack:${TORCH_OVERLAY}/site-packages:${BASE}/site-packages" \
            "cu130-base-standalone:${BASE}/site-packages"; do
  label="${cell%%:*}"
  path="${cell#*:}"
  PYTHONPATH="${path}:${ROOT}" python -m unittest -v \
    tests.unit.test_b6k_admission_probe \
    > "${STAGE}/unittest-${label}.stdout.txt" \
    2> "${STAGE}/unittest-${label}.stderr.txt"
  PYTHONPATH="${path}:${ROOT}" python - \
    "${COBOLT_CKPT}" "${COBOLT_SHA256}" "${STAGE}/probe3-${label}.json" <<'PY'
from pathlib import Path
import json
import sys

from scripts.b6k_admission_probe import probe_checkpoint_load
import torch

record = probe_checkpoint_load(Path(sys.argv[1]), sys.argv[2])
record["torch_version"] = torch.__version__
record["torch_file"] = str(Path(torch.__file__).resolve())
Path(sys.argv[3]).write_text(
    json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
if record["status"] != "pass":
    raise SystemExit(f"probe 3 did not pass: {record['status']}")
print(json.dumps({k: record[k] for k in ("status", "tensor_count", "state_digest")}))
PY
done

# The strict-restore digest must not depend on which torch performed the load.
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

stage = Path(sys.argv[1])
records = {p.name: json.loads(p.read_text()) for p in sorted(stage.glob("probe3-*.json"))}
digests = {name: value["state_digest"] for name, value in records.items()}
if len(set(digests.values())) != 1:
    raise SystemExit(f"state digest differs across runtimes: {digests}")
(stage / "cross_runtime_agreement.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-b6k-probe3-cross-runtime-v1",
            "runtimes": {n: v["torch_version"] for n, v in records.items()},
            "state_digest": next(iter(digests.values())),
            "agree": True,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
print(json.dumps(digests, indent=2))
PY

sha256sum \
  "${ROOT}/scripts/b6k_admission_probe.py" \
  "${ROOT}/tests/unit/test_b6k_admission_probe.py" \
  "${ROOT}/slurm/preflight_b6k_probe3_checkpoint_cpu.sbatch" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]).resolve(strict=True),
    {
        "artifact_class": "b6k_preflight_probe3_checkpoint",
        "gpu_submitted": False,
        "dispatcher_files_touched": False,
        "outcomes_read": False,
        "project_data_read": False,
        "metrics_calculated": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'preflight artifact: %s\n' "${TARGET}"
