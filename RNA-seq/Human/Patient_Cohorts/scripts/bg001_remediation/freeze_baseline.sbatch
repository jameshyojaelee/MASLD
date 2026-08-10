#!/usr/bin/env bash
#SBATCH --job-name=edgeR
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
: "${BG001_RUN_ROOT:?BG001_RUN_ROOT is required}"
[[ -f "$BG001_RUN_ROOT/.bg001_candidate_root" ]] || { echo "Missing candidate sentinel" >&2; exit 64; }
SCRIPT_DIR="$BG001_RUN_ROOT/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
mkdir -p "$BG001_RUN_ROOT/logs"
exec > >(tee "$BG001_RUN_ROOT/logs/freeze_baseline.out") 2> >(tee "$BG001_RUN_ROOT/logs/freeze_baseline.err" >&2)
eval "$(micromamba shell hook --shell=bash)"
set +u
micromamba activate rnaseq
set -u
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$BG001_RUN_ROOT" --require-bam-frozen
python3 "$SCRIPT_DIR/verify_analysis_contract.py" --run-root "$BG001_RUN_ROOT"
python3 - "$BG001_RUN_ROOT" "$PROJECT_ROOT" <<'PY'
import json
import pathlib
import sys

run_root = pathlib.Path(sys.argv[1]).resolve(strict=True)
project_root = pathlib.Path(sys.argv[2]).resolve(strict=True)
contract_root = pathlib.Path(
    json.loads((run_root / "contract/run_contract.json").read_text())["project_root"]
).resolve(strict=True)
if project_root != contract_root:
    raise SystemExit(
        "MASLD_PROJECT_ROOT differs from frozen run-contract project_root: "
        f"{project_root} != {contract_root}"
    )
PY
python3 "$SCRIPT_DIR/analysis_runtime_contract.py" create --run-root "$BG001_RUN_ROOT"
python3 "$SCRIPT_DIR/analysis_runtime_contract.py" verify --run-root "$BG001_RUN_ROOT"
R_SCRIPT=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["Rscript"]["path"])' \
  "$BG001_RUN_ROOT/contract/analysis_runtime_contract.json")
"$R_SCRIPT" --vanilla "$SCRIPT_DIR/freeze_baseline.R" "$PROJECT_ROOT" "$BG001_RUN_ROOT"
python3 "$SCRIPT_DIR/analysis_runtime_contract.py" verify --run-root "$BG001_RUN_ROOT"
python3 "$SCRIPT_DIR/verify_run_contract.py" --run-root "$BG001_RUN_ROOT" --require-bam-frozen
python3 "$SCRIPT_DIR/seal_baseline.py" --run-root "$BG001_RUN_ROOT"
python3 "$SCRIPT_DIR/verify_analysis_contract.py" \
  --run-root "$BG001_RUN_ROOT" --require-baseline-frozen
