#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-matrix-preflight
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ROSTER="${ROOT}/config/campaigns/sequence_native_diagonal_input_matrix_20260824.tsv"
BUILDER="${ROOT}/slurm/build_chrombpnet_hepatocyte_donor0_genomic0_inputs.sbatch"
PSEUDOBULK="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
SPLITS="${ROOT}/executions/sequence-split-contract-v2-21063827"
TARGET="${ROOT}/executions/sequence-native-diagonal-input-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
bash -n "${BUILDER}"
printf '%s  %s\n' 3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b \
  "${PSEUDOBULK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLITS}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROSTER}" "${PSEUDOBULK}" "${SPLITS}" \
  > "${STAGE}/validation/preflight.json" <<'PY'
from pathlib import Path
import csv
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

roster_path, pseudobulk, splits = map(Path, sys.argv[1:])
verify_frozen_tree(pseudobulk)
verify_frozen_tree(splits)
with roster_path.open(encoding="utf-8", newline="") as handle:
    roster = list(csv.DictReader(handle, delimiter="\t"))
with (pseudobulk / "training_pseudobulk_manifest.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    pseudobulk_rows = list(csv.DictReader(handle, delimiter="\t"))
with (splits / "crossed_outer_splits.tsv").open(encoding="utf-8", newline="") as handle:
    split_rows = {row["split_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
if len(roster) != 24 or len({(r["lineage_id"], r["outer_fold"]) for r in roster}) != 24:
    raise SystemExit("diagonal input roster differs")
for row in roster:
    outer = int(row["outer_fold"])
    if outer not in range(5) or row["split_id"] != f"donor{outer}_genomic{outer}":
        raise SystemExit("diagonal roster identity differs")
    matches = [
        item for item in pseudobulk_rows
        if int(item["donor_test_fold"]) == outer
        and item["lineage_id"] == row["lineage_id"]
    ]
    if len(matches) != 1 or int(matches[0]["donors"]) < 1:
        raise SystemExit("pseudobulk roster binding differs")
    split = split_rows.get(row["split_id"])
    if split is None or int(split["genomic_valid_fold"]) != (outer + 1) % 5:
        raise SystemExit("diagonal split binding differs")
print(json.dumps({
    "schema_version": "masld-bench-sequence-native-diagonal-input-preflight-v1",
    "status": "pass", "campaign_rows": len(roster),
    "lineages": sorted({row["lineage_id"] for row in roster}),
    "outer_folds": sorted({int(row["outer_fold"]) for row in roster}),
    "pseudobulk_frozen": True, "splits_frozen": True,
    "completed_hepatocyte_donor0_genomic0_excluded": True,
}, sort_keys=True, indent=2))
PY
sha256sum "${ROSTER}" "${BUILDER}" \
  "${ROOT}/slurm/test_sequence_native_diagonal_input_matrix.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/preflight.json").read_text())
if receipt["status"] != "pass" or receipt["campaign_rows"] != 24:
    raise SystemExit("diagonal preflight receipt differs")
freeze_tree(root, {"artifact_class": "sequence_native_diagonal_input_preflight", **receipt})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
