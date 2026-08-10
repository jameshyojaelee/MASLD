#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

python3 "$SCRIPT_DIR/tests/test_genetics_preflight.py"
if command -v Rscript >/dev/null 2>&1; then
  Rscript "$SCRIPT_DIR/tests/test_power_fixture.R"
else
  micromamba run -n rnaseq Rscript "$SCRIPT_DIR/tests/test_power_fixture.R"
fi
