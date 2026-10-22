"""lncRNA inset geometry. Left: forest of unadjusted and adjusted odds ratios (log10 x axis; ticks 0.5 and 1.0; the
gray reference line sits at 1.0). Right top: share up/down per biotype (y ticks 0%, 50%, 100%). Right bottom: matched
direction odds ratio on a log10 axis (ticks 0.1, 0.3, 1). Compared with our recomputation: unadjusted 1.119
(1.007-1.245), adjusted 0.429 (0.377-0.487; job 9834), shares 83.0% / 47.9%, and the matched OR reported in RESULTS
0.132 (0.090-0.193), which we reproduce only in direction (0.145-0.193 under our own matching). Usage: <pdf> <tag>"""
import pdfplumber, numpy as np, pandas as pd, sys
pdf,tag=sys.argv[1],sys.argv[2]; p=pdfplumber.open(pdf).pages[0]
W=pd.DataFrame(p.extract_words()); W['cx']=(W.x0+W.x1)/2; W['cy']=(W.top+W.bottom)/2
col=lambda c:tuple(round(v,4) for v in (c if isinstance(c,(list,tuple)) else [c]*3))
ticks=[l for l in p.lines if abs(l['x0']-l['x1'])<0.01 and 1<(l['bottom']-l['top'])<2]
def tl(l):
    c=W[(abs(W.cx-l['x0'])<6)&(W.top>l['bottom'])&(W.top<l['bottom']+8)&W.text.str.fullmatch(r'[0-9.]+')]; return float(c.text.iloc[0]) if len(c) else None
xt=[(l['x0'],tl(l)) for l in ticks if tl(l) is not None]
left=[t for t in xt if t[0]<145]; right=[t for t in xt if t[0]>190]
print("left ticks",left,"right ticks",right)
fl=np.polyfit([t[0] for t in left],np.log10([t[1] for t in left]),1); fr=np.polyfit([t[0] for t in right],np.log10([t[1] for t in right]),1)
FL=lambda x:10**np.polyval(fl,x); FR=lambda x:10**np.polyval(fr,x)
out=[]
# forest rows by label position
rowlab={'Unadjusted':W[W.text=='Unadjusted'].cy.iloc[0],'Adjusted':W[W.text=='Adjusted'].cy.iloc[0]}
vlines=[l for l in p.lines if abs(l['x0']-l['x1'])<0.01 and (l['bottom']-l['top'])>3 and (l['bottom']-l['top'])<10 and col(l['stroking_color'])==(0,0,0)]
hlines=[l for l in p.lines if abs(l['top']-l['bottom'])<0.01 and (l['x1']-l['x0'])>1 and l['top']<140 and col(l['stroking_color'])==(0,0,0) and l['x1']-l['x0']<60]
marks=[q for q in p.rects if 2<(q['x1']-q['x0'])<6 and 2<(q['bottom']-q['top'])<6 and q['x1']<145]+[c for c in p.curves if (c['x1']-c['x0'])<6]
exp={'Unadjusted':(1.119,1.007,1.245),'Adjusted':(0.429,0.377,0.487)}
for k,yl in rowlab.items():
    ys=[l for l in hlines if l['x1']<145 and abs(l['top']-yl)<12]; m=[q for q in marks if q['x1']<145 and abs((q['top']+q['bottom'])/2-yl)<12]
    if m: out.append(dict(mark='forest_point',key=k,drawn=FL((m[0]['x0']+m[0]['x1'])/2),expected=exp[k][0]))
    if ys: out.append(dict(mark='forest_ci_low',key=k,drawn=FL(ys[0]['x0']),expected=exp[k][1])); out.append(dict(mark='forest_ci_high',key=k,drawn=FL(ys[0]['x1']),expected=exp[k][2]))
# bars
yt=[l for l in p.lines if abs(l['top']-l['bottom'])<0.01 and 1<(l['x1']-l['x0'])<2 and l['x0']>150]
def yl_(l): c=W[(abs(W.cy-l['top'])<2.5)&(W.x1<l['x0']+0.5)&(W.x1>l['x0']-15)]; return float(c.text.iloc[0].replace('%','')) if len(c) else None
ytv=[(l['top'],yl_(l)) for l in yt if yl_(l) is not None]; fy=np.polyfit([t[0] for t in ytv],[t[1] for t in ytv],1); print("bar y ticks",ytv)
MAG=(0.7882,0.149,0.3686)
cats={w:W[(W.text==w)&(W.top>100)&(W.x0>150)].cx.iloc[0] for w in ['Protein-coding','lncRNA'] if ((W.text==w)&(W.top>100)&(W.x0>150)).any()}
shares={'Protein-coding':83.0,'lncRNA':47.9}
for q in p.rects:
    if col(q['non_stroking_color'])==MAG and q['x0']>150 and (q['bottom']-q['top'])>10 and q['top']>30:
        c=min(cats,key=lambda k:abs(cats[k]-(q['x0']+q['x1'])/2)); out.append(dict(mark='bar_up_share',key=c,drawn=np.polyval(fy,q['top'])-np.polyval(fy,q['bottom']),expected=shares[c]))
# matched
mm=[c for c in p.curves if c['x0']>190 and (c['x1']-c['x0'])<6]
if mm: out.append(dict(mark='matched_point',key='matched',drawn=FR((mm[0]['x0']+mm[0]['x1'])/2),expected=0.132))
mv=[l for l in vlines if l['x0']>190 and l['top']>110]; mh=[l for l in hlines if l['x0']>190 and l['top']>110]
print("matched: whisker caps drawn at",[round(FR(l['x0']),4) for l in mv],"| horizontal CI segments drawn:",len(mh),"| axis left edge value",round(FR(min(l['x0'] for l in p.lines if l['x0']>190 and abs(l['top']-l['bottom'])<0.01 and (l['x1']-l['x0'])>50)),4))
O=pd.DataFrame(out); O['off']=[abs(np.log10(d)-np.log10(e))>0.005*abs(np.log10(1.3)-np.log10(0.075)) if m.startswith('forest') or m.startswith('matched') else abs(d-e)>0.5 for m,d,e in zip(O.mark,O.drawn,O.expected)]
print(O.round(4).to_string()); print("marks checked",len(O),"off",int(O.off.sum()))
O.to_csv(f"out/lncrna_geometry_{tag}.tsv",sep="\t",index=False)
