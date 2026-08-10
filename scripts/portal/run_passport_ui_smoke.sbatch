#!/bin/bash
#SBATCH --job-name=portal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/chromium_smoke_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/chromium_smoke_%j.err

set -euo pipefail

module load Mesa/23.1.9-GCCcore-13.2.0

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
UI_SMOKE_DIR="$PROJECT_ROOT/RNA-seq/results/evidence_passports/candidates/program-context-v2-candidate-2026-08-07/fixtures/ui_chromium_smoke"
UI_CHROME="${PASSPORT_CHROMIUM:-/gpfs/commons/home/jameslee/.cache/puppeteer/chrome/linux-150.0.7871.24/chrome-linux64/chrome}"

mkdir -p "$UI_SMOKE_DIR"
"$UI_CHROME" \
  --headless=new \
  --no-sandbox \
  --disable-gpu \
  --hide-scrollbars \
  --window-size=800,600 \
  --screenshot="$UI_SMOKE_DIR/smoke.png" \
  "data:text/html,%3Chtml%3E%3Cbody%20style%3D%27font-family%3Asans-serif%27%3E%3Ch1%3EPassport%20UI%20renderer%3C%2Fh1%3E%3Cp%3ELocal-only%20smoke%20test%3C%2Fp%3E%3C%2Fbody%3E%3C%2Fhtml%3E"

"$UI_CHROME" --version
file "$UI_SMOKE_DIR/smoke.png"
stat -c "%s bytes" "$UI_SMOKE_DIR/smoke.png"
