#!/usr/bin/env bash
set -euo pipefail

MODE=${1:-}
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SRC=${ROOT}/GWAS/finemapping/src/seqfunc_v2/mpra
V2=${ROOT}/GWAS/finemapping/results/seqfunc/mpra_benchmark/v2

case "${MODE}" in
  source_env)
    source_job=$(sbatch --parsable --job-name=MPRA "${SRC}/02_source_reproduction.sbatch")
    env_job=$(sbatch --parsable --job-name=MPRAnalyze "${SRC}/00_build_env.sbatch")
    gate_job=$(sbatch --parsable --job-name=MPRA --dependency="afterany:${source_job}:${env_job}" "${SRC}/04_gate.sbatch")
    echo "source_job=${source_job}"
    echo "environment_job=${env_job}"
    echo "gate_job=${gate_job}"
    ;;
  smoke)
    python3 - "${V2}/gate_verdict.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]))
if not (x["source_reproduction_pass"] and x["mpranalyze_environment_complete"]):
    raise SystemExit("Source/environment gates do not permit MPRAnalyze smoke submission")
PY
    smoke_job=$(sbatch --parsable --job-name=MPRAnalyze "${SRC}/03_mpranalyze_smoke.sbatch")
    gate_job=$(sbatch --parsable --job-name=MPRA --dependency="afterany:${smoke_job}" "${SRC}/04_gate.sbatch")
    echo "smoke_job=${smoke_job}"
    echo "gate_job=${gate_job}"
    ;;
  full)
    python3 - "${V2}/gate_verdict.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]))
if not x["mpranalyze_full_submission_allowed"]:
    raise SystemExit("Source/environment/smoke gates do not permit full MPRAnalyze submission")
PY
    full_job=$(sbatch --parsable --job-name=MPRAnalyze "${SRC}/03_mpranalyze.sbatch")
    gate_job=$(sbatch --parsable --job-name=MPRA --dependency="afterany:${full_job}" "${SRC}/04_gate.sbatch")
    echo "full_job=${full_job}"
    echo "gate_job=${gate_job}"
    ;;
  *)
    echo "Usage: $0 {source_env|smoke|full}" >&2
    exit 2
    ;;
esac
