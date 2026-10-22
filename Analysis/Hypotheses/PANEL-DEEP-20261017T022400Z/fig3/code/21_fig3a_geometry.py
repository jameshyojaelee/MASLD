"""Fig 3A geometry: bar segments (disease, control) per cohort row from the x-axis tick transform, printed totals, and
the 15 coverage tiles (Sex, NAS, Fibrosis) decoded through the legend swatches. Compared with out/fig3a_recomputed.tsv.
Usage: <pdf> <tag>"""
import pdfplumber, numpy as np, pandas as pd, sys
pdf,tag=sys.argv[1],sys.argv[2]; p=pdfplumber.open(pdf).pages[0]
W=pd.DataFrame(p.extract_words()); R=pd.read_csv("out/fig3a_recomputed.tsv",sep="\t").set_index('dataset')
col=lambda q:tuple(round(v,4) for v in (q['non_stroking_color'] if isinstance(q['non_stroking_color'],(list,tuple)) else [q['non_stroking_color']]*3))
rows=W[W.text.str.startswith('GSE')&(W.x0<45)].copy(); rows['cy']=(rows.top+rows.bottom)/2
# x axis of bar plot: tick lines below the right plot, labels directly beneath
xt=[l for l in p.lines if abs(l['x0']-l['x1'])<0.01 and (l['bottom']-l['top'])<2 and l['x0']>145]
def lab(x): c=W[(abs((W.x0+W.x1)/2-x)<4)&(W.top>xt[0]['bottom'])&(W.top<xt[0]['bottom']+8)]; return float(c.text.iloc[0])
tv=[(l['x0'],lab(l['x0'])) for l in xt]; a=np.polyfit([t[0] for t in tv],[t[1] for t in tv],1); print("x ticks",tv)
DIS=(0.7882,0.149,0.3686); CTL=(0.6196,0.6196,0.6196)
bars=[q for q in p.rects if q['x0']>145 and (q['bottom']-q['top'])>5 and col(q) in (DIS,CTL) and q['top']<100]
out=[]
for _,r in rows.iterrows():
    segs=[q for q in bars if q['top']<=r.cy<=q['bottom']]
    for q in segs:
        arm='disease' if col(q)==DIS else 'control'; v=(q['x1']-q['x0'])*a[0]
        out.append(dict(mark='bar',key=f"{r.text}:{arm}",drawn=v,expected=R.loc[r.text,arm],start=np.polyval(a,q['x0'])))
    num=W[(W.x0>145)&(abs((W.top+W.bottom)/2-r.cy)<2)&W.text.str.fullmatch(r'\d+')]
    if len(num): out.append(dict(mark='label',key=f"{r.text}:total",drawn=float(num.text.iloc[0]),expected=R.loc[r.text,'n']))
# legend swatches -> state
leg={}
for w in ['Absent','Partial','Binned','Available']:
    ww=W[W.text==w].iloc[0]; sw=[q for q in p.rects if q['x1']<=ww.x0+0.5 and q['x1']>ww.x0-8 and abs((q['top']+q['bottom'])/2-(ww.top+ww.bottom)/2)<3 and col(q)!=(1,1,1)]
    leg[col(sw[-1])]=w
print("legend",leg)
def expstate(d,var):
    n=R.loc[d,'n']
    if var=='Sex': k=R.loc[d,'sex_inferred_nonNA']
    elif var=='NAS': k=R.loc[d,'nas_nonNA']
    else:
        if d=='GSE213621': return 'Binned'
        k=R.loc[d,'fib_nonNA']
    return 'Available' if k==n else ('Absent' if k==0 else 'Partial')
tiles=[q for q in p.rects if q['x1']<145 and q['x0']>40 and (q['bottom']-q['top'])>10 and col(q) in leg]
xs=sorted(set(round(q['x0'],1) for q in tiles)); vars_=['Sex','NAS','Fibrosis']
for q in tiles:
    cy=(q['top']+q['bottom'])/2; d=rows.iloc[(rows.cy-cy).abs().argmin()].text; v=vars_[xs.index(round(q['x0'],1))]
    out.append(dict(mark='tile',key=f"{d}:{v}",drawn=leg[col(q)],expected=expstate(d,v)))
O=pd.DataFrame(out)
O['off']=[(str(a)!=str(b)) if m=='tile' else abs(float(a)-float(b))>max(0.5,0.005*(tv[-1][1]-tv[0][1])) for m,a,b in zip(O.mark,O.drawn,O.expected)]
print(O.to_string()); print("marks checked",len(O),"off",int(O.off.sum()))
print("axis title words:",W[(W.top>xt[0]['bottom']+8)&(W.x0>145)&(W.top<125)].text.tolist())
O.to_csv(f"out/fig3a_geometry_{tag}.tsv",sep="\t",index=False)
