"""Fig 3C geometry: recover study-bar, overlap-bar and ribbon widths (memberships) from the vector paths and compare with
the recomputed cohort x class membership matrix. Usage: python 14_fig3c_geometry.py <pdf> <tag>"""
import pdfplumber, pandas as pd, numpy as np, sys
pdf,tag=sys.argv[1],sys.argv[2]
p=pdfplumber.open(pdf).pages[0]
M=pd.read_csv("out/fig3c_recomputed_cohort_by_class.tsv",sep="\t",index_col=0)
rects=[r for r in p.rects if r['fill']]
left=sorted([r for r in rects if r['x1']-r['x0']<12 and r['x0']<40],key=lambda r:r['top'])
right=sorted([r for r in rects if r['x1']-r['x0']<12 and r['x0']>100],key=lambda r:r['top'])
tot=sum(r['bottom']-r['top'] for r in left); k=M.loc['All','All']/tot
print("scale memberships/pt",round(k,3))
# study order by label positions: words 'GSE' followed by number; leader lines map label to bar
W=p.extract_words()
nums=[w for w in W if w['text'].isdigit() and len(w['text'])==6]
lines=p.lines
def bar_for_label(w):
    y=(w['top']+w['bottom'])/2
    # leader line starting near label
    for l in lines:
        if l['x0']<25 and abs(l['top']-y)<6:
            yend=l['bottom']; return min(left,key=lambda r:0 if r['top']<=yend<=r['bottom'] else 1e9)
    return min(left,key=lambda r:0 if r['top']<=y<=r['bottom'] else 1e9)
rows=[]
for w in nums:
    b=bar_for_label(w); ds='GSE'+w['text']
    rows.append(dict(mark='study_bar',key=ds,drawn=(b['bottom']-b['top'])*k,expected=M.loc[ds,'All'],color=str(b['non_stroking_color'])))
classes=['unique','same2','same3','same4','same5','mixed']
for r,c in zip(right,classes):
    rows.append(dict(mark='overlap_bar',key=c,drawn=(r['bottom']-r['top'])*k,expected=M.loc['All',c],color=str(r['non_stroking_color'])))
# ribbons: curves + full-width rects; left span -> study, right span -> class
rib=[]
for c in p.curves:
    xs=[q[0] for q in c['pts']]; L=[q[1] for q in c['pts'] if abs(q[0]-min(xs))<0.01]; R=[q[1] for q in c['pts'] if abs(q[0]-max(xs))<0.01]
    rib.append((min(L),max(L),min(R),max(R),c['non_stroking_color']))
for r in rects:
    if r['x1']-r['x0']>50: rib.append((r['top'],r['bottom'],r['top'],r['bottom'],r['non_stroking_color']))
study_of={}
for w in nums: study_of[id(bar_for_label(w))]='GSE'+w['text']
for (l0,l1,r0,r1,col) in rib:
    lb=[b for b in left if b['top']-0.01<=(l0+l1)/2<=b['bottom']+0.01][0]; rb=[b for b in right if b['top']-0.01<=(r0+r1)/2<=b['bottom']+0.01][0]
    ds=study_of[id(lb)]; cl=classes[right.index(rb)]
    rows.append(dict(mark='ribbon',key=f"{ds}->{cl}",drawn=(l1-l0)*k,drawn_right=(r1-r0)*k,expected=M.loc[ds,cl],color=str(col)))
O=pd.DataFrame(rows); O['off']=(O.drawn-O.expected).abs()>max(0.5,0.005*M.loc['All','All'])
print(O.round(1).to_string()); print("marks checked",len(O),"off",O.off.sum())
O.to_csv(f"out/fig3c_geometry_{tag}.tsv",sep="\t",index=False)
# printed class labels
txt=' '.join(w['text'] for w in W); print(txt)
