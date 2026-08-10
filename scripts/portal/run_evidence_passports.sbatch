#!/bin/bash
#SBATCH --job-name=portal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/passport_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/passport_%j.err

set -euo pipefail

module load Mesa/23.1.9-GCCcore-13.2.0

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
# `sbatch` executes a private spool copy, so BASH_SOURCE cannot identify the
# submitted file's original siblings.  Normal runs use the repository tools;
# replay from a REL01 archive must set PASSPORT_TOOL_ROOT to the frozen
# `environments/producers/portal` directory and therefore cannot silently fall
# back to mutable producers.
TOOL_ROOT="${PASSPORT_TOOL_ROOT:-$PROJECT_ROOT/scripts/portal}"
GENERATOR="$TOOL_ROOT/generate_evidence_passports.py"
VALIDATOR="$TOOL_ROOT/validate_evidence_passports.py"
BUILDER="$TOOL_ROOT/build_evidence_passport_bundle.py"
FIXTURE_TESTER="$TOOL_ROOT/test_evidence_passport_bundle.py"
CANDIDATE_ROOT="${PASSPORT_CANDIDATE_ROOT:-$PROJECT_ROOT/RNA-seq/results/evidence_passports/candidates/program-context-v2-candidate-2026-08-07}"
PASSPORT_MODE="${PASSPORT_MODE:-production}"
PASSPORT_SELECTION="${PASSPORT_SELECTION:-$CANDIDATE_ROOT/passport_input_selection.tsv}"
PASSPORT_FINALIZE="${PASSPORT_FINALIZE:-false}"

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PASSPORT_MESA_MODULE="${PASSPORT_MESA_MODULE:-Mesa/23.1.9-GCCcore-13.2.0}"

if [[ "$PASSPORT_MODE" == "fixtures" ]]; then
  micromamba run -n spatial python -B "$GENERATOR" \
    --mode fixtures --candidate-root "$CANDIDATE_ROOT"
  micromamba run -n spatial python -B "$VALIDATOR" \
    --self-test --candidate-root "$CANDIDATE_ROOT"
  micromamba run -n spatial python -B "$FIXTURE_TESTER"
elif [[ "$PASSPORT_MODE" == "production" ]]; then
  if [[ "$PASSPORT_FINALIZE" == "true" ]]; then
    # Seal the exact preterminal bytes that were manually reviewed.  Do not
    # rebuild HTML or screenshots here: the v2 attestation binds their hashes.
    if [[ -z "${PASSPORT_MANUAL_ACCEPTANCE:-}" ]]; then
      echo "PASSPORT_MANUAL_ACCEPTANCE is required when PASSPORT_FINALIZE=true" >&2
      exit 2
    fi
    micromamba run -n spatial python -B "$VALIDATOR" \
      --candidate-root "$CANDIDATE_ROOT" --selection "$PASSPORT_SELECTION" \
      --bundle "$CANDIDATE_ROOT" --finalize-pass06 \
      --manual-acceptance "$PASSPORT_MANUAL_ACCEPTANCE"
  else
    # Every build command fails closed until the signed coordinator selection,
    # terminal GEN closure, and real Plan13 READY seal match selected artifacts.
    micromamba run -n spatial python -B "$GENERATOR" \
      --mode production --candidate-root "$CANDIDATE_ROOT" \
      --selection "$PASSPORT_SELECTION"
    micromamba run -n spatial python -B "$BUILDER" \
      --selection "$PASSPORT_SELECTION" --output-root "$CANDIDATE_ROOT"
    micromamba run -n spatial python -B "$VALIDATOR" \
      --candidate-root "$CANDIDATE_ROOT" --selection "$PASSPORT_SELECTION" \
      --bundle "$CANDIDATE_ROOT"
  fi
elif [[ "$PASSPORT_MODE" == "verify-sealed" ]]; then
  # This branch is deliberately read-only and accepts a relocated bundle root.
  micromamba run -n spatial python -B "$VALIDATOR" \
    --bundle "$CANDIDATE_ROOT" --verify-sealed
else
  echo "PASSPORT_MODE must be production, fixtures, or verify-sealed" >&2
  exit 2
fi
