#!/usr/bin/env bash
#SBATCH --job-name=data-check-244
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=executions/logs/data-check-244-%j.out
#SBATCH --error=executions/logs/data-check-244-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
FAILED_CAMPAIGN=${ROOT}/executions/model-train-243-21092120-0.failed/campaign
TASK_ID=gse274114_hiseq_healthy_vs_hbv
VIEW_SHA256=(
  5603f451d196970a4ce79aea89d2a3ff5a36770bbec613220df25266e99445c1
  e070358394bd226868669acf433ba89cc7bbc3bc8cf641af0a9d4e499d358a77
  c2bb539456877457fd158605a3c654805d87fb96d26091b6ae4bb7d0d380ff66
  1d29b7a33a445eeca74af943b2e7f6da25b6447bd056fc249e51c0c9ee283c41
  f2682d04aa1716235106200e10f2ddfa31a982e9fb21ee1203d27326fe2d3265
)
TARGET=${ROOT}/executions/data-check-244-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE274114 convergence audit.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
for fold in 0 1 2 3 4; do
  printf '%s  %s\n' "${VIEW_SHA256[${fold}]}" \
    "${FAILED_CAMPAIGN}/trainer_views/fold-${fold}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
"${PYTHON}" - <<'PY'
from pathlib import Path

path = Path("scripts/diagnose_gse274114_convergence.py")
compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
bash -n slurm/diagnose_gse274114_convergence_cpu.sbatch

"${PYTHON}" scripts/diagnose_gse274114_convergence.py \
  --task-id "${TASK_ID}" \
  --view-root "${FAILED_CAMPAIGN}/trainer_views/fold-0" \
  --view-root "${FAILED_CAMPAIGN}/trainer_views/fold-1" \
  --view-root "${FAILED_CAMPAIGN}/trainer_views/fold-2" \
  --view-root "${FAILED_CAMPAIGN}/trainer_views/fold-3" \
  --view-root "${FAILED_CAMPAIGN}/trainer_views/fold-4" \
  --view-sha256 "${VIEW_SHA256[0]}" \
  --view-sha256 "${VIEW_SHA256[1]}" \
  --view-sha256 "${VIEW_SHA256[2]}" \
  --view-sha256 "${VIEW_SHA256[3]}" \
  --view-sha256 "${VIEW_SHA256[4]}" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/audit.stdout.json" \
  2> "${STAGE}/audit.stderr.txt"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
manifest = verify_frozen_tree(root / "audit")
receipt = json.loads((root / "audit/receipt.json").read_text(encoding="utf-8"))
with (root / "audit/selected_candidate_summary.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if (
    (root / "audit.stderr.txt").read_text(encoding="utf-8")
    or manifest["metadata"].get("artifact_class")
    != "gse274114_training_only_convergence_audit"
    or receipt.get("status") != "passed_training_only_convergence_audit"
    or receipt.get("diagnostic_runs") != 50
    or receipt.get("query_molecular_values_opened") is not False
    or receipt.get("query_labels_opened") is not False
    or receipt.get("predictions_generated") is not False
    or receipt.get("scoring_or_evaluation_run") is not False
    or len(rows) != 50
):
    raise SystemExit("training-only convergence audit contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse274114_training_only_convergence_audit_execution",
        "status": receipt["status"],
        "diagnostic_runs": 50,
        "query_labels_opened": False,
        "predictions_generated": False,
        "scoring_or_evaluation_run": False,
    },
)
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/diagnose_gse274114_convergence.py \
  scripts/gse274114_baseline_firewall.py \
  slurm/diagnose_gse274114_convergence_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/audit/ARTIFACTS.json"
