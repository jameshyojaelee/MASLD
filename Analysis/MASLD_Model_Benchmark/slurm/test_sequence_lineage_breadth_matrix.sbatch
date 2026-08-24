#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-breadth-gate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ROSTER="${ROOT}/config/campaigns/sequence_native_lineage_breadth_screen_matrix_20260824.tsv"
PRIOR_GATE="${ROOT}/executions/sequence-hepatocyte-screen-gate-21066927"
TARGET="${ROOT}/executions/sequence-lineage-breadth-gate-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 73f76b29518353bacab55e74e7be7b37eea448d54497f01028cde8fa7c4e609f \
  "${PRIOR_GATE}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROSTER}" "${PRIOR_GATE}" \
  > "${STAGE}/validation/gate.json" <<'PY'
from pathlib import Path
import csv
import json
import sys
from collections import Counter
from masld_bench.artifacts import verify_frozen_tree

roster_path, prior_gate = map(Path, sys.argv[1:])
verify_frozen_tree(prior_gate)
with roster_path.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
models = {"sequence_cnn_control", "sequence_transformer_control", "bpnet", "chrombpnet"}
lineages = {"cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell"}
keys = {(r["model_id"], r["lineage_id"], int(r["outer_fold"]), int(r["seed"])) for r in rows}
if len(rows) != 120 or len(keys) != 120:
    raise SystemExit("lineage breadth matrix cardinality differs")
for row in rows:
    fold = int(row["outer_fold"])
    seed = int(row["seed"])
    if (
        row["model_id"] not in models
        or row["lineage_id"] not in lineages
        or fold not in range(5)
        or seed not in {20260824, 20260825, 20260826}
        or row["split_id"] != f"donor{fold}_genomic{fold}"
    ):
        raise SystemExit("lineage breadth row identity differs")
seed_counts = Counter(int(row["seed"]) for row in rows)
existing = [row for row in rows if row["disposition"].startswith("completed_existing_")]
submitted = [row for row in rows if row["disposition"] == "submit"]
if seed_counts != {20260824: 100, 20260825: 10, 20260826: 10}:
    raise SystemExit("lineage breadth seed allocation differs")
if len(existing) != 4 or len(submitted) != 116:
    raise SystemExit("lineage breadth disposition differs")
print(json.dumps({
    "schema_version": "masld-bench-sequence-lineage-breadth-gate-v1",
    "status": "pass", "matrix_rows": len(rows), "submission_rows": len(submitted),
    "existing_rows": len(existing), "seed_counts": dict(sorted(seed_counts.items())),
    "lineages": sorted(lineages), "models": sorted(models), "folds": list(range(5)),
}, sort_keys=True, indent=2))
PY
sha256sum "${ROSTER}" "${ROOT}/slurm/test_sequence_lineage_breadth_matrix.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/gate.json").read_text())
freeze_tree(root, {"artifact_class": "sequence_lineage_breadth_gate", **receipt})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
