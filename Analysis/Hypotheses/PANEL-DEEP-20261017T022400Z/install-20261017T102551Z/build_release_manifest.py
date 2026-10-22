"""Release manifest: one release ID over every paper output file, with sha256. Usage: build_release_manifest.py RELEASE_ID OUTDIR"""
import os, sys, hashlib, csv, glob
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'; os.chdir(ROOT)
RID, OUT = sys.argv[1], sys.argv[2]; os.mkdir(OUT)
CAT='Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda: f.read(1<<20), b''): h.update(b)
    return h.hexdigest()
rows=[]
def add(p, role):
    rows.append({'release_id':RID,'role':role,'path':p,'bytes':os.path.getsize(p),'sha256':sha(p)})
for dp,_,fs in os.walk('figures/main'):
    if '/_archive' in dp or '/candidates' in dp: continue
    for f in sorted(fs): add(os.path.join(dp,f), 'figure_pdf' if f.endswith('.pdf') else 'figure_source_or_index')
for dp,_,fs in os.walk('docs/manuscript'):
    if dp.startswith('docs/manuscript/reviews'): continue
    for f in sorted(fs): add(os.path.join(dp,f),'manuscript_text_or_ledger')
for f in ['README.md','PAPER.md','STATUS.md','RESULTS.md','ROADMAP.md']: add('docs/'+f,'authority')
for f in sorted(glob.glob(CAT+'/*')): add(f,'gene_catalog_r9d')
with open(OUT+'/release_manifest.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
from collections import Counter
print(RID, len(rows), dict(Counter(r['role'] for r in rows)))
