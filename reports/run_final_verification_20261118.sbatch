#!/usr/bin/env bash
#SBATCH --job-name=final_verify_261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/reports/final_verify_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/reports/final_verify_20261118_%j.err

set -euo pipefail

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${ROOT}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PY_SP="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
PY_RN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
PY_AGA="/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python"
RS_RN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

echo "=== [1/9] Resource Scope Check (--strict-latex) ==="
python3 scripts/manuscript/validate_resource_scope.py --strict-latex

echo "=== [2/9] Suite 1: Portal, Manuscript, Figure & Core Analysis Unit Tests (spatial) ==="
"${PY_SP}" -B -m pytest -p no:cacheprovider -q \
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

echo "=== [3/9] Suite 2: alt_contig_recount Unit Tests (rnaseq w/ pysam) ==="
"${PY_RN}" -B -m pytest -p no:cacheprovider -q \
  scripts/analysis/alt_contig_recount/test_annotation_correspondence.py \
  scripts/analysis/alt_contig_recount/test_count_execution_evidence.py \
  scripts/analysis/alt_contig_recount/test_native_library_subsample.py \
  scripts/analysis/alt_contig_recount/test_verified_fragments.py

echo "=== [4/9] Suite 3: AlphaGenome Atlas Unit Tests (alphagenome_atlas) ==="
"${PY_AGA}" -B -m unittest discover -s scripts/analysis/alphagenome_atlas/tests -p "test_*.py"

echo "=== [5/9] Suite 4: ATAC Regulatory Methods Unit Tests (spatial) ==="
"${PY_SP}" -B -m pytest -p no:cacheprovider -q Analysis/ATAC/Regulatory_Methods/tests/

echo "=== [6/9] Suite 5: R Statistical & Identity Unit Tests (rnaseq) ==="
for rf in \
  Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/test_identity_guard.R \
  scripts/manuscript/resource_f_five/test_qualified_bulk_inputs.R \
  scripts/figures/test_pi_stage_inputs.R
do
  echo "Running R test: ${rf}"
  "${RS_RN}" --vanilla "${rf}"
done

echo "=== [7/9] Model .1 Release Verification (3 test_release.sh scripts) ==="
bash Analysis/MASLD_Model_Benchmark/release/masld-severity-v1.2.1/test_release.sh
bash Analysis/MASLD_Model_Benchmark/release/masld-panel-v1.1/test_release.sh
bash Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2.1/test_release.sh

echo "=== [8/9] N22 Spot Check ==="
"${PY_SP}" -B - <<'PYEOF'
import sys, inspect, numpy as np
sys.path.insert(0, "Analysis/MASLD_Model_Benchmark/executions/standardized-metric-panel-20260902T130923Z")
from panel_metrics import midrank
assert np.allclose(midrank([3.0, 1.0, 1.0, 2.0]), [4.0, 1.5, 1.5, 3.0])
try:
    midrank([1.0, np.nan, 2.0])
    raise AssertionError("midrank accepted NaN")
except ValueError:
    print("PASS: panel_metrics.midrank raises ValueError on NaN")

sys.path.insert(0, "scripts/analysis/masld_two_models")
import chromatin_development as cd
assert inspect.signature(cd.fit_nuisance_ridge).parameters["allow_zero_ridge"].default is False
print("PASS: chromatin_development.fit_nuisance_ridge defaults allow_zero_ridge=False")
PYEOF

echo "=== [9/9] Programmatic Remediation Verification (N01-N23) ==="
python3 reports/verify_remediation_20261118.py

echo "=== FINAL_VERIFICATION_OVERALL_PASS ==="
