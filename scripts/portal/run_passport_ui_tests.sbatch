#!/bin/bash
#SBATCH --job-name=portal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/passport_ui_test_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/passport_ui_test_%j.err

set -euo pipefail

module load Mesa/23.1.9-GCCcore-13.2.0

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PASSPORT_MESA_MODULE="Mesa/23.1.9-GCCcore-13.2.0"
export PASSPORT_CHROMIUM="${PASSPORT_CHROMIUM:-/gpfs/commons/home/jameslee/.cache/puppeteer/chrome/linux-150.0.7871.24/chrome-linux64/chrome}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

cd "$PROJECT_ROOT"

micromamba run -n spatial python -B -m py_compile \
  scripts/portal/build_evidence_passport_ui.py \
  scripts/portal/render_evidence_passport_ui.py \
  scripts/portal/build_evidence_passport_bundle.py \
  scripts/portal/generate_evidence_passports.py \
  scripts/portal/prepare_passport_attestations.py \
  scripts/portal/validate_evidence_passports.py \
  scripts/portal/test_evidence_passport_bundle.py \
  scripts/manuscript/program_context_v2/rel01_snapshot_candidate.py \
  scripts/manuscript/program_context_v2/publish_plan60_workstream_handoff.py \
  scripts/manuscript/program_context_v2/coordinator/coordinator_contract.py \
  scripts/manuscript/program_context_v2/coordinator/tests/test_coordinator_contract.py

micromamba run -n spatial ruff check --select E9,F63,F7,F82 \
  scripts/portal/build_evidence_passport_ui.py \
  scripts/portal/render_evidence_passport_ui.py \
  scripts/portal/build_evidence_passport_bundle.py \
  scripts/portal/generate_evidence_passports.py \
  scripts/portal/prepare_passport_attestations.py \
  scripts/portal/validate_evidence_passports.py \
  scripts/portal/test_evidence_passport_bundle.py

micromamba run -n spatial python -B scripts/portal/generate_evidence_passports.py \
  --mode fixtures \
  --candidate-root RNA-seq/results/evidence_passports/candidates/program-context-v2-candidate-2026-08-07
micromamba run -n spatial python -B scripts/portal/validate_evidence_passports.py \
  --self-test \
  --candidate-root RNA-seq/results/evidence_passports/candidates/program-context-v2-candidate-2026-08-07
micromamba run -n spatial python -B scripts/portal/test_evidence_passport_bundle.py

micromamba run -n spatial bash scripts/manuscript/program_context_v2/tests/run_tests.sh
micromamba run -n spatial bash scripts/manuscript/program_context_v2/coordinator/tests/run_tests.sh

git diff --check -- \
  scripts/portal/build_evidence_passport_ui.py \
  scripts/portal/render_evidence_passport_ui.py \
  scripts/portal/build_evidence_passport_bundle.py \
  scripts/portal/generate_evidence_passports.py \
  scripts/portal/prepare_passport_attestations.py \
  scripts/portal/validate_evidence_passports.py \
  scripts/portal/test_evidence_passport_bundle.py \
  scripts/portal/run_evidence_passports.sbatch \
  scripts/portal/run_passport_ui_tests.sbatch \
  scripts/manuscript/program_context_v2
