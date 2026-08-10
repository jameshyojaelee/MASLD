#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

python3 "$SCRIPT_DIR/tests/test_release_scaffold.py"
