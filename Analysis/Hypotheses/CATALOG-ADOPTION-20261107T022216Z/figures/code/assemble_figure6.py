"""Current Figure6 review: same centered vertical order, native scale, enough height for current panels."""
from pathlib import Path
import csv, hashlib,json,subprocess
from collections import Counter
from pypdf import PdfReader,PdfWriter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
B=ROOT/'Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures'
font_manager.fontManager.addfont('/gpfs/commons/home/jameslee/.local/share/fonts/Helvetica.ttf')
plt.rcParams.update({'font.family':'Helvetica','font.size':6,'pdf.fonttype':42})
panels=[(p.name[4].upper(),p) for p in sorted((B/'panels').glob('fig6*.pdf'))]
readers=[(letter,p,PdfReader(p)) for letter,p in panels]
W=468; gap=18; top=16; H=top+sum(float(r.pages[0].mediabox.height) for _,_,r in readers)+gap*(len(readers)-1)+8
w=PdfWriter();page=w.add_blank_page(width=W,height=H);rows=[]
fig=plt.figure(figsize=(W/72,H/72),facecolor='none');ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,W),ylim=(0,H));ax.axis('off')
for letter,path,reader in readers:
 p=reader.pages[0];pw,ph=float(p.mediabox.width),float(p.mediabox.height);x=(W-pw)/2;y=H-top-ph
 assert x>=0 and y>=0
 page.merge_translated_page(p,x,y);ax.text(x,y+ph+3,letter,fontsize=6,va='bottom')
 rows.append(dict(callout=letter,source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),width_pt=pw,height_pt=ph,x_pt=x,top_pt=top,scale=1));top+=ph+gap
ov=B/'review_assembly/Figure6_letters.pdf';fig.savefig(ov,transparent=True);plt.close(fig);page.merge_page(PdfReader(ov).pages[0]);page.compress_content_streams()
out=B/'review_assembly/Figure6.pdf';assert not out.exists()
with out.open('wb') as f:w.write(f)
def text(p):return subprocess.check_output(['pdftotext','-layout',str(p),'-'],text=True)
a=text(out);expected=Counter()
for _,p in panels:expected+=Counter(''.join(text(p).split()))
assert not(expected-Counter(''.join(a.split())))
with (B/'review_assembly/Figure6_panel_sources.tsv').open('w') as f:
 d=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');d.writeheader();d.writerows(rows)
(B/'review_assembly/Figure6_checks.json').write_text(json.dumps(dict(pages=1,width_pt=W,height_pt=H,scale=1,panel_characters_missing=0,layout='centered vertical native panels, same order'),indent=2)+'\n')
print(out)
