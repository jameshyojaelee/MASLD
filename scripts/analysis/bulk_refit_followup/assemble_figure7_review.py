"""Review proposal from six indexed Figure7 PDFs, translation only."""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path

repo=Path.cwd(); out=Path(os.environ['FIG34_OUT']); (out/'proofs').mkdir()
spec=importlib.util.spec_from_file_location('native_assembler',repo/'figures/candidates/fig3-fig4-assembled-774-20261013T074505Z/producers/assemble_fig34.py')
assembler=importlib.util.module_from_spec(spec);spec.loader.exec_module(assembler)
index=repo/'figures/main/fig7_model_framework/CANDIDATE_PANEL_INDEX.tsv'
rows=list(csv.DictReader(index.open(),delimiter='\t'));assert len(rows)==6
paths=[repo/r['candidate_relative_path'] for r in rows]
sizes=[assembler.size(p) for p in paths]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
margin=18.;gap=24.
column_widths=[max(sizes[i][0] for i in range(c,6,2)) for c in range(2)]
row_heights=[max(sizes[i][1] for i in range(r*2,r*2+2)) for r in range(3)]
assembler.PAGE_W=2*margin+sum(column_widths)+gap
height=2*margin+sum(row_heights)+2*gap
slots=[];basis={};top=height-margin
for i,(path,(w,h)) in enumerate(zip(paths,sizes)):
    r,c=divmod(i,2)
    if c==0 and r: top-=row_heights[r-1]+gap
    x=margin+(column_widths[0]+gap if c else 0)
    letter=rows[i]['callout'][-1]
    slots.append(dict(letter=letter,x=x,y=top-h,w=w,h=h,path=path,letter_x=x,letter_y=top+12,note=None))
    basis[('Figure7_review_proposal',letter)]='Current indexed source PDF; native-size translation-only review proposal; source claims not reinterpreted'
manifest,check=assembler.assemble('Figure7','Figure7_review_proposal',height,slots,basis)
assert hashes=={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
assert check['panel_characters_missing_from_composite_layout_text']==0
assert not check['type3_fonts']
# A review assembly preserves native panel fonts; inherited unembedded fonts
# are reported rather than modifying source artwork.
source_fonts=set()
for path in paths:
    source_fonts.update(k[0] for k,v in assembler.fonts_of(assembler.PdfReader(path)).items() if not v)
assert set(check['not_embedded']) <= source_fonts
with (out/'panel_manifest.tsv').open('w') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(manifest[0]),delimiter='\t');writer.writeheader();writer.writerows(manifest)
(out/'source_identity.json').write_text(json.dumps({'source_sha256_before_equals_after':hashes,'assembly_sha256':hashlib.sha256((out/'Figure7_review_proposal.pdf').read_bytes()).hexdigest(),'proposal_only':True},indent=2)+'\n')
print(json.dumps(check,indent=2))
