#!/usr/bin/env bash
set -uo pipefail

if [[ "$(hostname)" == *login* ]]; then
  echo "ERROR: Refusing to run unit tests on login node $(hostname)" >&2
  exit 1
fi

REPO="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
REPRO="$REPO/reports/codebase-deep-review-20261118/repro"
LOG="$REPRO/unit_test_suites.log"
TMPRUN="$REPRO/tmp_pytest"
mkdir -p "$TMPRUN"
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$TMPRUN"

MM="/gpfs/commons/home/jameslee/.local/bin/micromamba"
cd "$REPO"

{
  echo "=== Compute Node Unit Test Verification ==="
  echo "Hostname: $(hostname)"
  echo "Slurm JobId: ${SLURM_JOB_ID:-NONE}"
  echo "UTC Start: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo ""

  echo "--- Suite 1: Resource Scope Validator (rapids_singlecell, python3) ---"
  t0=$(date +%s.%N)
  "$MM" run -n rapids_singlecell python3 -B scripts/manuscript/validate_resource_scope.py
  rc1=$?
  t1=$(date +%s.%N)
  echo "ExitCode=$rc1 Elapsed=$(python3 -c "print(round($t1 - $t0, 2))")s"
  echo ""

  echo "--- Suite 2: Portal, Manuscript, Figure & Core Analysis Python Unit Tests (spatial, pytest) ---"
  t0=$(date +%s.%N)
  "$MM" run -n spatial python3 -B -m pytest -p no:cacheprovider --continue-on-collection-errors -q \
    scripts/portal/test_evidence_passport_bundle.py \
    scripts/portal/test_gene_catalog_observability.py \
    scripts/portal/test_gene_catalog_v2.py \
    scripts/portal/test_spatial_resource_semantics.py \
    scripts/analysis/fig34_corrections/tests/test_resource_identity.py \
    scripts/analysis/masld_two_models/test_scientific_invariants.py \
    scripts/figures/test_fig2_locus_coordinates.py \
    scripts/figures/tests/test_fig2_cluster_inference.py \
    scripts/manuscript/aging_inspired_resource_figures/tests/test_gnmt_state.py \
    scripts/manuscript/noncoding_resource/tests/test_annotate_program_lncrna_content.py \
    scripts/manuscript/noncoding_resource/tests/test_audit_strandedness.py \
    scripts/manuscript/noncoding_resource/tests/test_gene_identity.py \
    scripts/manuscript/noncoding_resource/tests/test_noncoding_genetics.py \
    scripts/manuscript/noncoding_resource/tests/test_validate_program_without_lncrna.py \
    scripts/manuscript/program_context_v2/coordinator/tests/test_coordinator_contract.py \
    scripts/manuscript/program_context_v2/tests/test_release_scaffold.py \
    Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/test_identity_aggregation.py \
    Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/xenium/code/test_selection_replay.py
  rc2=$?
  t1=$(date +%s.%N)
  echo "ExitCode=$rc2 Elapsed=$(python3 -c "print(round($t1 - $t0, 2))")s"
  echo ""

  echo "--- Suite 3: AlphaGenome Atlas Python Unit Tests (alphagenome2, unittest) ---"
  t0=$(date +%s.%N)
  "$MM" run -n alphagenome2 python3 -B -m unittest discover -s scripts/analysis/alphagenome_atlas/tests -p "test_*.py" -v
  rc3=$?
  t1=$(date +%s.%N)
  echo "ExitCode=$rc3 Elapsed=$(python3 -c "print(round($t1 - $t0, 2))")s"
  echo ""

  echo "--- Suite 4: ATAC Regulatory Methods Python Unit Tests (spatial, pytest) ---"
  t0=$(date +%s.%N)
  "$MM" run -n spatial python3 -B -m pytest -p no:cacheprovider --continue-on-collection-errors -q \
    Analysis/ATAC/Regulatory_Methods/tests/
  rc4=$?
  t1=$(date +%s.%N)
  echo "ExitCode=$rc4 Elapsed=$(python3 -c "print(round($t1 - $t0, 2))")s"
  echo ""

  echo "--- Suite 5: R Statistical & Identity Unit Tests (rnaseq, Rscript) ---"
  t0=$(date +%s.%N)
  r_pass=0
  r_fail=0
  for rf in \
    Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/test_identity_guard.R \
    scripts/manuscript/resource_f_five/test_qualified_bulk_inputs.R \
    scripts/figures/test_pi_stage_inputs.R
  do
    echo "Running R test: $rf"
    if "$MM" run -n rnaseq Rscript --vanilla "$rf"; then
      r_pass=$((r_pass + 1))
    else
      r_fail=$((r_fail + 1))
    fi
  done
  t1=$(date +%s.%N)
  echo "R_Tests_Passed=$r_pass R_Tests_Failed=$r_fail Elapsed=$(python3 -c "print(round($t1 - $t0, 2))")s"
  echo ""
  echo "UTC End: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} 2>&1 | tee "$LOG"

rm -rf "$TMPRUN"
find "$REPO/Analysis/MASLD_Model_Benchmark/release" -name "*.pyc" -newer "$REPRO/../snapshot_utc.txt" -delete 2>/dev/null || true
