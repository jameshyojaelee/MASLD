#!/usr/bin/env python3
"""Translation-only review proposal from indexed Figure 2 A-I PDFs."""
import csv,hashlib,json,os
from pathlib import Path
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from PyPDF2 import PageObject,PdfReader,PdfWriter,Transformation
ROOT=Path(os.environ['MASLD_PROJECT_ROOT']); OUT=Path(os.environ['FIG2_PROPOSAL_ROOT']); OUT.mkdir(parents=True,exist_ok=False)
BASE=ROOT/'figures/main/fig2_genetics/panels'
with (BASE/'Figure2_reindex_manifest.tsv').open() as fh:rows=list(csv.DictReader(fh,delimiter='\t'))
selected={x['panel']:x for x in rows if x['panel'] in 'ABCDEFGHI'}
assert set(selected)==set('ABCDEFGHI')
# Native geometry: three-panel rows, two locus examples, shared locus key.
W,H=7.2,9.2
positions={'A':(.185,.25),'B':(3.915,.25),'C':(5.915,.25),'D':(.185,2.90),'E':(2.235,2.90),'F':(4.435,2.90),'G':(.185,5.30),'H':(3.695,5.30),'I':(.185,8.70)}
plt.rcParams.update({'pdf.fonttype':42,'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],'font.size':6})
fig=plt.figure(figsize=(W,H));ax=fig.add_axes([0,0,1,1]);ax.set_axis_off()
for name,(x,top) in positions.items():ax.text(x/W,1-(top-.015)/H,name,ha='left',va='bottom',fontsize=6,fontweight='normal',transform=ax.transAxes)
label=OUT/'panel_labels.pdf';fig.savefig(label,facecolor='white');plt.close(fig)
page=PdfReader(str(label)).pages[0];records=[];bounds=[]
for name,(x,top) in positions.items():
 p=BASE/selected[name]['canonical_pdf'];sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha==selected[name]['sha256'],name
 reader=PdfReader(str(p));assert len(reader.pages)==1;src=reader.pages[0]
 w,h=float(src.mediabox.width)/72,float(src.mediabox.height)/72
 assert x+w<=W+1e-8 and top+h<=H+1e-8
 for old,ox,ot,ow,oh in bounds:assert x+w<=ox or ox+ow<=x or top+h<=ot or ot+oh<=top,(name,old)
 bounds.append((name,x,top,w,h));ctm=Transformation().translate(tx=x*72,ty=(H-top-h)*72).ctm
 assert ctm[:4]==(1,0,0,1)
 # Same transformation-aware merge used by the existing Figure 5 utility;
 # transform the clipping rectangle and content together, never rescale.
 page._merge_page(src,lambda content,src=src,ctm=ctm:PageObject._add_transformation_matrix(content,src.pdf,ctm),ctm,False)
 records.append(dict(panel=name,source=str(p),source_sha256=sha,width_in=w,height_in=h,x_in=x,top_in=top,transformation=list(ctm)))
writer=PdfWriter();writer.add_page(page)
with (OUT/'Figure2_A-I_native_review_PROPOSAL.pdf').open('wb') as fh:writer.write(fh)
with (OUT/'panel_identity_geometry.json').open('w') as fh:json.dump(dict(page_width_in=W,page_height_in=H,panels=records,scientific_state='PROPOSAL: indexed panels unchanged; final journal sizing open'),fh,indent=2)
companions=[str(BASE/x['canonical_pdf']) for x in rows if x['panel']=='J']+[str(BASE/'Fig2_companion_noncoding_annotation_gap.pdf'),str(BASE/'Fig2_companion_eqtl_target_gene_biotype.pdf')]
(OUT/'unassembled_companion_references.txt').write_text('\n'.join(companions)+'\n')
(OUT/'CHECKED.txt').write_text('All nine source hashes match current index. All single-page native boxes fit and do not overlap. Every placement has identity linear transform and translation only. Individual source PDFs unchanged. Proposal only; journal sizing remains open.\n')
print(OUT)
