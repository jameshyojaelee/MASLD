"""Fig 3G geometry: parse the content stream (matplotlib marker XObjects placed by cumulative 'cm' translations, colour
set by 'g'/'rg'), recover every point's (log2FC, PP.H4) from the tick transform, and compare with the recomputed table per
class via nearest-neighbour matching both ways. Usage: <pdf> <tag>"""
import re, sys, numpy as np, pandas as pd, pdfplumber
from pypdf import PdfReader
from scipy.spatial import cKDTree
pdf,tag=sys.argv[1],sys.argv[2]
pg=PdfReader(pdf).pages[0]; H=float(pg.mediabox.height)
c=pg.get_contents().get_data().decode('latin1')
c=re.sub(r'BT.*?ET',' ',c,flags=re.S)
toks=re.findall(r'/[A-Za-z0-9]+|[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?|[A-Za-z*\']+',c)
stack=[]; ctm=np.eye(3); fill=(0,0,0); st=[]; pts=[]
i=0; ops=[]
nums=[]
for t in toks:
    if re.fullmatch(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',t): nums.append(float(t)); continue
    if t.startswith('/'): nums.append(t); continue
    op=t
    if op=='q': st.append((ctm.copy(),fill))
    elif op=='Q': ctm,fill=st.pop()
    elif op=='cm':
        a,b,cc,d,e,f=nums[-6:]; M=np.array([[a,b,0],[cc,d,0],[e,f,1]]); ctm=M@ctm
    elif op=='rg' and len(nums)>=3: fill=tuple(nums[-3:])
    elif op=='g' and nums: fill=(nums[-1],)*3
    elif op=='Do': pts.append((ctm[2,0],ctm[2,1],tuple(round(v,4) for v in fill),nums[-1]))
    nums=[]
P=pd.DataFrame(pts,columns=['x','y','fill','xobj']); P['top']=H-P.y
# axes from pdfplumber ticks
pp=pdfplumber.open(pdf).pages[0]
xt=sorted([l for l in pp.lines if abs(l['x1']-l['x0'])<0.01 and l['bottom']-l['top']<3],key=lambda l:l['x0'])
yt=sorted([l for l in pp.lines if abs(l['bottom']-l['top'])<0.01 and l['x1']-l['x0']<3],key=lambda l:l['top'])
W=pd.DataFrame(pp.extract_words())
def lab(x=None,y=None):
    if x is not None: r=W[(abs((W.x0+W.x1)/2-x)<4)&(W.top>xt[0]['bottom'])&(W.top<xt[0]['bottom']+10)]
    else: r=W[(abs((W.top+W.bottom)/2-y)<3)&(W.x1<yt[0]['x0']+1)]
    return float(r.text.iloc[0].replace('−','-'))
xv=[(l['x0'],lab(x=l['x0'])) for l in xt]; yv=[(l['top'],lab(y=l['top'])) for l in yt]
ax=np.polyfit(*zip(*xv),1); ay=np.polyfit(*zip(*yv),1)
print("x ticks",xv,"y ticks",yv)
# marker xobject centre offset: compute from xobject bbox
xo=pg['/Resources']['/XObject']
off={}
for k in xo:
    bb=[float(v) for v in xo[k].get_object()['/BBox']]; off[k]=((bb[0]+bb[2])/2,(bb[1]+bb[3])/2)
P['cx']=P.x+[off[k][0] for k in P.xobj]; P['cy']=P.top-[off[k][1] for k in P.xobj]
P['lfc']=np.polyval(ax,P.cx); P['pp4']=np.polyval(ay,P.cy)
print(P.groupby('fill').size())
# legend markers are outside plotting area (below axis)
inplot=P[P.cy<88.5]; leg=P[P.cy>=88.5]
cmap={}
R=pd.read_csv("out/fig3g_recomputed.tsv",sep="\t")
for f,g in inplot.groupby('fill'):
    best=None
    for cl,h in R.groupby('cls'):
        if len(h)==len(g): best=cl
    cmap[f]=best
print("colour->class by count:",cmap)
rows=[]; rng_x=np.ptp(np.polyval(ax,[27,188.6])); rng_y=1.0
for f,g in inplot.groupby('fill'):
    h=R[R.cls==cmap[f]]
    A=np.c_[g.lfc/rng_x,g.pp4/rng_y]; B=np.c_[h.logFC/rng_x,h.coloc_best_susie_pp4/rng_y]
    d1,_=cKDTree(B).query(A); d2,_=cKDTree(A).query(B)
    tol=0.005
    rows.append(dict(fill=f,cls=cmap[f],n_drawn=len(g),n_expected=len(h),max_nn_drawn_to_exp=d1.max(),max_nn_exp_to_drawn=d2.max(),off_drawn=(d1>tol).sum(),off_exp=(d2>tol).sum()))
O=pd.DataFrame(rows); print(O.to_string())
print("x data range drawn",inplot.lfc.min(),inplot.lfc.max(),"expected",R.logFC.min(),R.logFC.max())
P.to_csv(f"out/fig3g_drawn_points_{tag}.tsv",sep="\t",index=False); O.to_csv(f"out/fig3g_geometry_{tag}.tsv",sep="\t",index=False)
