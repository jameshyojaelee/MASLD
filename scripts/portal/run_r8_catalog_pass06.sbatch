#!/bin/bash
#SBATCH --job-name=r8-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/r8_review_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/r8_review_%j.err
set -euo pipefail
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
TASK_REVIEW_ROOT="$TASK_ROOT/figures/candidates/catalog-r8-reviewed-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
TASK_SOURCE="$TASK_ROOT/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z"
TASK_CLONE="$TASK_REVIEW_ROOT/catalog_reviewed_r8"
TASK_PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
export MASLD_PROJECT_ROOT="$TASK_ROOT"
export PASSPORT_CANDIDATE_ROOT="$TASK_CLONE"
export PASSPORT_ANALYSIS_RELEASE_ID=gene-catalog-corrected-r8-20261009T142027Z
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$TASK_ROOT"
mkdir "$TASK_REVIEW_ROOT"
printf '%s\n' "$TASK_REVIEW_ROOT" > "$TASK_ROOT/logs/r8_review_${SLURM_JOB_ID}.run_root.txt"
cp -a "$TASK_SOURCE" "$TASK_CLONE"
TASK_SOURCE="$TASK_SOURCE" TASK_CLONE="$TASK_CLONE" TASK_REVIEW_ROOT="$TASK_REVIEW_ROOT" "$TASK_PYTHON" -B - <<'PY'
import hashlib,json,os
from pathlib import Path
source=Path(os.environ['TASK_SOURCE']);clone=Path(os.environ['TASK_CLONE']);out=Path(os.environ['TASK_REVIEW_ROOT'])
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
before={str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file()}
assert before=={str(p.relative_to(clone)):sha(p) for p in clone.rglob('*') if p.is_file()}
(out/'exact_clone_checksums.json').write_text(json.dumps(before,indent=2)+'\n')
(out/'actual_scientific_review.json').write_text(json.dumps({'reviewer':'Codex remaining_catalog scientific review','reviewer_is_user':False,'original_images_inspected':['overview','THRB','HKDC1','GLP1R','MTARC1'],'observed_classes':{'THRB':'genetically_anchored','HKDC1':'concordant','GLP1R':'untested','MTARC1':'unresolved'},'visible_checks':['gene/program/assay grain boundary','source dependence','testability distinction','next experiment and falsifier'],'limitation':'Agent scientific review, not a forged user attestation. Historical Myojin class assays remain stale and cannot qualify r8 classes. Bulk inferential units are conservative genetic-identity clusters, not verified independent donors; 855 held annotations remain untestable. No external publication.'},indent=2)+'\n')
PY
"$TASK_PYTHON" -B scripts/portal/prepare_passport_attestations.py manual \
 --selection "$TASK_CLONE/passport_input_selection.tsv" --bundle "$TASK_CLONE" \
 --output "$TASK_REVIEW_ROOT/passport_manual_acceptance.tsv" \
 --review-id "r8-scientific-review-20261014-${SLURM_JOB_ID}" \
 --reviewer 'Codex remaining_catalog scientific review' \
 --reviewed-at-utc "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
 --decision-register-id 'user-r8-selection-20261010-and-complete-resource-20261014' \
 --confirm-call-states --confirm-provenance-states --confirm-visible-boundary \
 --confirm-visible-testability --confirm-visible-source-dependence --confirm-visible-falsifier
"$TASK_PYTHON" -B scripts/portal/validate_evidence_passports.py \
 --candidate-root "$TASK_CLONE" --selection "$TASK_CLONE/passport_input_selection.tsv" \
 --bundle "$TASK_CLONE" --finalize-pass06 --manual-acceptance "$TASK_REVIEW_ROOT/passport_manual_acceptance.tsv"
"$TASK_PYTHON" -B scripts/portal/validate_evidence_passports.py --bundle "$TASK_CLONE" --verify-sealed
TASK_SOURCE="$TASK_SOURCE" TASK_CLONE="$TASK_CLONE" TASK_REVIEW_ROOT="$TASK_REVIEW_ROOT" "$TASK_PYTHON" -B - <<'PY'
import hashlib,json,os
from pathlib import Path
source=Path(os.environ['TASK_SOURCE']);clone=Path(os.environ['TASK_CLONE']);out=Path(os.environ['TASK_REVIEW_ROOT'])
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
before=json.loads((out/'exact_clone_checksums.json').read_text())
assert before=={str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file()}
scientific=['passport_gene_index.parquet','passport_evidence_long.parquet','passport_coverage_long.parquet','passport_program_index.parquet','passport_program_membership.parquet','passport_program_context.parquet','passport_assay_status.tsv','passport_input_selection.tsv','passport_input_selection.signature.json','portal_candidate/index.html','portal_candidate/review/review_manifest.tsv']
assert all(sha(clone/p)==before[p] for p in scientific)
(out/'review_completion.json').write_text(json.dumps({'status':'PASS','reviewed_clone':str(clone),'historical_source_modified':False,'scientific_tables_unchanged':scientific,'manual_reviewer':'Codex remaining_catalog scientific review','external_publication':False,'input_source_files_unchanged':len(before)},indent=2)+'\n')
PY
