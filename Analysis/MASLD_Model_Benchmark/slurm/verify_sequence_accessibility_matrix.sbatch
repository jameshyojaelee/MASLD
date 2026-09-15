#!/usr/bin/env bash
#SBATCH --job-name=seq-eval-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MATRIX="${ROOT}/executions/sequence-accessibility-valid-matrix-21070421"
MATRIX_SHA=7960cb07a58f7d43d66615da87bc6347d94d6e488d763414678049480fa38fe3
TARGET="${ROOT}/executions/sequence-accessibility-valid-matrix-audit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite sequence matrix audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${MATRIX_SHA}" "${MATRIX}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${MATRIX}" "${STAGE}" "${MATRIX_SHA}" <<'PY'
from collections import Counter
import csv
import json
from pathlib import Path
import re
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

matrix, output = Path(sys.argv[1]).resolve(strict=True), Path(sys.argv[2]).resolve(strict=True)
verify_frozen_tree(matrix)
summary = json.loads((matrix / "matrix_summary.json").read_text(encoding="utf-8"))
if (
    summary["status"] != "pass"
    or summary["evaluation_groups"] != 7
    or summary["evaluation_role"] != "valid"
    or summary["folds_or_seeds_used_as_biological_replicates"]
    or summary["cross_fold_candidate_ranking_calculated"]
    or summary["incomplete_coverage_used_for_comparison"]
    or summary["test_outcomes_read"]
    or summary["sealed_data_read"]
    or summary["champion_claim_allowed"]
):
    raise SystemExit("matrix summary safety contract differs")
candidate_rows = 0
unit_rows = 0
role_counts = Counter()
for evaluation in summary["evaluations"]:
    result_path = matrix / evaluation["result"]
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if (
        result["evaluation_role"] != "valid"
        or result["biological_unit"] != "donor"
        or result["genomic_unit"] != "held_contig_block"
        or result["test_outcomes_read"]
        or result["champion_claim_allowed"]
        or result["external_evaluation"]
        or result["seeds_used_as_biological_replicates"]
    ):
        raise SystemExit(f"evaluation contract differs: {result_path}")
    candidate_rows += len(result["prediction_artifacts"])
    metrics_path = result_path.parent / "donor_block_metrics.tsv"
    with metrics_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"evaluation has no donor/block units: {metrics_path}")
    for row in rows:
        role_counts[row["role"]] += 1
        if (
            row["role"] != "valid"
            or re.fullmatch(r"[0-9a-f]{64}", row["donor_hash"]) is None
            or re.fullmatch(r"[0-9a-f]{64}", row["block_hash"]) is None
        ):
            raise SystemExit(f"metric unit identity differs: {metrics_path}")
    unit_rows += len(rows)
audit = {
    "schema_version": "masld-bench-sequence-accessibility-valid-matrix-audit-v1",
    "status": "pass",
    "matrix_artifacts_sha256": sys.argv[3],
    "evaluation_groups": len(summary["evaluations"]),
    "candidate_evaluation_rows": candidate_rows,
    "donor_block_metric_rows": unit_rows,
    "metric_role_counts": dict(sorted(role_counts.items())),
    "frozen_tree_verified": True,
    "raw_donor_ids_exported": False,
    "test_outcomes_read": False,
    "sealed_data_read": False,
    "champion_claim_allowed": False,
}
(output / "audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
freeze_tree(
    output,
    {
        "artifact_class": "sequence_accessibility_valid_matrix_audit",
        "matrix_artifacts_sha256": "7960cb07a58f7d43d66615da87bc6347d94d6e488d763414678049480fa38fe3",
        "evaluation_groups": audit["evaluation_groups"],
        "frozen_tree_verified": True,
        "test_outcomes_read": False,
        "sealed_data_read": False,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Sequence evaluation matrix audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
