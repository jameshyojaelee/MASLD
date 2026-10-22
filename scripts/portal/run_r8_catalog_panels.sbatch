#!/bin/bash
#SBATCH --job-name=r8-fig6
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/r8_fig6_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/r8_fig6_%j.err
set -euo pipefail
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
TASK_SOURCE="$TASK_ROOT/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects"
TASK_OUTPUT="$TASK_ROOT/figures/candidates/catalog-r8-panels-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
TASK_PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
mkdir "$TASK_OUTPUT"
printf '%s\n' "$TASK_OUTPUT" > "$TASK_ROOT/logs/r8_fig6_${SLURM_JOB_ID}.run_root.txt"
cp "$TASK_SOURCE/sidecar_completion.json" "$TASK_SOURCE/worked_entries.tsv" "$TASK_OUTPUT/"
cp "$TASK_ROOT/scripts/portal/render_resource_catalog_review.py" "$TASK_OUTPUT/render_resource_catalog_review.py"
cp "$TASK_ROOT/scripts/portal/run_r8_catalog_panels.sbatch" "$TASK_OUTPUT/run_r8_catalog_panels.sbatch"
"$TASK_PYTHON" -B "$TASK_OUTPUT/render_resource_catalog_review.py" --sidecar "$TASK_OUTPUT" --hide-status-lines --class-header 'Gene evidence class'
for TASK_PDF in "$TASK_OUTPUT"/panels/*.pdf; do
 pdfinfo "$TASK_PDF"
 pdftotext -layout "$TASK_PDF" "${TASK_PDF%.pdf}.txt"
 gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "$TASK_PDF"
done
TASK_ROOT="$TASK_ROOT" TASK_SOURCE="$TASK_SOURCE" TASK_OUTPUT="$TASK_OUTPUT" "$TASK_PYTHON" -B - <<'PY'
import hashlib,json,os
from pathlib import Path
source=Path(os.environ['TASK_SOURCE']);out=Path(os.environ['TASK_OUTPUT'])
for f in ['sidecar_completion.json','worked_entries.tsv']:
 assert (out/f).read_bytes()==(source/f).read_bytes()
pdfs=list((out/'panels').glob('*.pdf'));assert len(pdfs)==3
assert all(b'/Subtype /Type3' not in p.read_bytes() and b'/Type /Font' in p.read_bytes() for p in pdfs)
c=(out/'panels/fig6c_worked_entries.txt').read_text()
assert all(x in c for x in ['GLP1R','HKDC1','MTARC1','THRB','Concordant','Untested','Unresolved','Genetically anchored'])
assert 'Candidate' not in c and 'bulk identity' not in c
metadata=json.loads((out/'sidecar_completion.json').read_text());assert metadata['object_rows']==28233 and metadata['bulk_analysis_unit_count']==774
(out/'panel_checks.json').write_text(json.dumps({'status':'PASS','source_sidecar':str(source),'source_primary_classes':metadata['primary_class_counts'],'pdfs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pdfs},'gene_entries':['GLP1R','HKDC1','MTARC1','THRB'],'output':'three individual PDFs; no figure assembly; no gene calls recomputed'},indent=2)+'\n')
PY
