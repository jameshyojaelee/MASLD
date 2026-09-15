#!/usr/bin/env bash
#SBATCH --job-name=model-check-353
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
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TOLERANCE="${ROOT}/config/b6k_numeric_tolerance_prespecification_20260825.json"
TARGET="${ROOT}/executions/b6k-admission-comparison-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"; FAILED="${TARGET}.failed"
preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve EXIT
[[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]] && { echo "refusing overwrite" >&2; exit 1; }
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' a26a6f3fe0d77ef12706038eaede32f222421154c0972e228b5238febeaa722f "${TOLERANCE}" | sha256sum --check --strict
module load python/3.12.3-GCCcore-13.3.0
cd "${ROOT}"
CU124_L40S="${ROOT}/executions/b6k-admission-cell-cu124-l40s-21106571/probe/numeric_fixture.pt"
CU128_L40S="${ROOT}/executions/b6k-admission-cell-cu128-l40s-21106657/probe/numeric_fixture.pt"
CU128_B6K="${ROOT}/executions/b6k-admission-cell-cu128-b6k-21106658/probe/numeric_fixture.pt"
CU130_B6K="${ROOT}/executions/b6k-admission-cell-cu130-b6k-21106570/probe/numeric_fixture.pt"

PYTHONPATH="${BASE}/site-packages:${ROOT}" python - \
  "${STAGE}" "${TOLERANCE}" "${CU124_L40S}" "${CU128_L40S}" "${CU128_B6K}" "${CU130_B6K}" <<'PY'
import json, sys
from pathlib import Path
from scripts.b6k_admission_probe import compare_fixtures

stage, tol_path = Path(sys.argv[1]), Path(sys.argv[2])
cu124_l40s, cu128_l40s, cu128_b6k, cu130_b6k = (Path(a) for a in sys.argv[3:7])
tol = json.loads(tol_path.read_text())
interchangeable = tol["thresholds"]["interchangeable_max_rel_difference"]
same_runtime_only = tol["thresholds"]["same_runtime_only_upper_bound"]

comparisons = {
    "architecture_effect__cu128_l40s_vs_cu128_b6k": (cu128_l40s, cu128_b6k),
    "version_effect__cu124_l40s_vs_cu128_l40s": (cu124_l40s, cu128_l40s),
    "combined__cu124_l40s_vs_cu128_b6k": (cu124_l40s, cu128_b6k),
    "b6k_runtime_effect__cu128_b6k_vs_cu130_b6k": (cu128_b6k, cu130_b6k),
}
results = {}
for name, (a, b) in comparisons.items():
    r = compare_fixtures(a, b)
    worst_rel = r.get("max_rel_difference", r.get("worst_relative"))
    worst_abs = r.get("max_abs_difference", r.get("worst_absolute"))
    if worst_rel is None:
        pt = r.get("per_tensor", {})
        vals = [v.get("max_rel_difference") for v in pt.values() if isinstance(v, dict) and "max_rel_difference" in v]
        worst_rel = max(vals) if vals else None
        avals = [v.get("max_abs_difference") for v in pt.values() if isinstance(v, dict) and "max_abs_difference" in v]
        worst_abs = max(avals) if avals else None
    verdict = "undetermined"
    if worst_rel is not None:
        if worst_rel <= interchangeable: verdict = "interchangeable"
        elif worst_rel <= same_runtime_only: verdict = "same_runtime_only_not_interchangeable"
        else: verdict = "materially_different"
    results[name] = {"max_rel_difference": worst_rel, "max_abs_difference": worst_abs, "verdict": verdict, "raw": r}
    (stage / f"{name}.json").write_text(json.dumps(r, indent=2, sort_keys=True) + "\n")

summary = {
    "schema_version": "masld-bench-b6k-admission-comparison-v1",
    "tolerance_prespecification_sha256": "a26a6f3fe0d77ef12706038eaede32f222421154c0972e228b5238febeaa722f",
    "thresholds": tol["thresholds"],
    "results": {k: {kk: vv for kk, vv in v.items() if kk != "raw"} for k, v in results.items()},
    "cluster_fact": (
        "L40S nodes ne1dg6-* run driver 550.107.02 / CUDA 12.4; B6K nodes ne1dg7-* run driver 590.48.01 / CUDA 13.1. "
        "torch 2.13.0+cu130 therefore cannot initialise on L40S (no_cuda_device_visible, verified job 21106572) and "
        "torch 2.6.0+cu124 cannot run on B6K (no_matching_sass, verified job 21106573). torch 2.8.0+cu128 is the ONLY "
        "runtime in the tree that executes on both pools."
    ),
    "outcomes_read": False, "metrics_calculated": False, "synthetic_inputs_only": True,
}
(stage / "comparison_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary["results"], indent=2, sort_keys=True))
PY

python3 -m pip freeze --all > "${STAGE}/pip_freeze.txt" 2>/dev/null || true
sha256sum "${ROOT}/slurm/compare_b6k_admission_cells_cpu.sbatch" "${TOLERANCE}" > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python - "${STAGE}" <<'PY'
import sys
from pathlib import Path
from masld_bench.artifacts import freeze_tree
freeze_tree(Path(sys.argv[1]), {"artifact_class": "b6k_admission_numeric_comparison",
  "synthetic_inputs_only": True, "biological_data_read": False,
  "development_outcomes_read": False, "sealed_outcomes_read": False,
  "benchmark_metrics_calculated": False})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
