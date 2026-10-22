"""Fig 3E geometry. Top: up (magenta) and down (blue) bars per stage from the y-tick transform (the lower tick is read as
negative because it sits below 0). Bottom: NMF activity points per program, with the program identified from the label
at the end of each line through its leader line. Compared with out/fig3e_recomputed_counts.tsv and _nmf.tsv.
Usage: <pdf> <tag>"""
import pdfplumber, numpy as np, pandas as pd, sys
pdf,tag=sys.argv[1],sys.argv[2]; p=pdfplumber.open(pdf).pages[0]
W=pd.DataFrame(p.extract_words(keep_blank_chars=False)); W['cy']=(W.top+W.bottom)/2; W['cx']=(W.x0+W.x1)/2
col=lambda c:tuple(round(v,4) for v in (c if isinstance(c,(list,tuple)) else [c]*3))
C=pd.read_csv("out/fig3e_recomputed_counts.tsv",sep="\t").set_index('contrast'); N=pd.read_csv("out/fig3e_recomputed_nmf.tsv",sep="\t")
yt=[l for l in p.lines if (l['x1']-l['x0'])<2 and abs(l['top']-l['bottom'])<0.01 and l['x0']<32]
def ylab(y): c=W[(abs(W.cy-y)<2.5)&(W.x1<30)&(W.upright)]; return c.text.iloc[0] if len(c) else None
ty=[(l['top'],ylab(l['top'])) for l in yt]; ty=[t for t in ty if t[1] is not None]
top=[t for t in ty if t[0]<50]; bot=[t for t in ty if t[0]>60]
y0=[t[0] for t in top if float(t[1])==0][0]
topv=[(y,float(s)*(1 if y<=y0 else -1)) for y,s in top]; signs_printed=[s for y,s in top]
a=np.polyfit([t[0] for t in topv],[t[1] for t in topv],1); b=np.polyfit([t[0] for t in bot],[float(t[1].replace('−','-')) for t in bot],1)
print("top ticks",top,"| bottom ticks",bot)
stx=W[W.text.str.fullmatch(r'F[0-4]')&(W.top<60)][['text','cx']].values.tolist()
stx2=W[W.text.str.fullmatch(r'F[0-4]')&(W.top>100)][['text','cx']].values.tolist()
MAG=(0.7882,0.149,0.3686); BLU=(0.0824,0.3961,0.7529)
out=[]
for q in p.rects:
    c=col(q['non_stroking_color']); 
    if c not in (MAG,BLU) or q['top']>60 or (q['x1']-q['x0'])>30: continue
    st=min(stx,key=lambda s:abs(s[1]-(q['x0']+q['x1'])/2))[0]; v=np.polyval(a,q['top'])-np.polyval(a,q['bottom'])
    key=f"{st}_vs_F0"; exp=C.loc[key,'up' if c==MAG else 'down']; out.append(dict(mark='bar',key=f"{st}:{'up' if c==MAG else 'down'}",drawn=abs(v),expected=exp,
        side='above0' if q['bottom']<=y0+0.01 else 'below0'))
rng_top=abs(topv[0][1]-topv[-1][1])
# NMF: labels P1..P5 and leader lines
labs=W[W.text.str.fullmatch(r'P[1-6]')]
leaders=[l for l in p.lines if l['x0']>150 and l['top']>60 and (l['x1']-l['x0'])>2]
pts=[c for c in p.curves if (c['x1']-c['x0'])<3 and c['top']>60]
def cyc(c): ys=[q[1] for q in c['pts']]; return (min(ys)+max(ys))/2
endF4=[(col(c['non_stroking_color']),cyc(c)) for c in pts if abs((c['x0']+c['x1'])/2-[s[1] for s in stx2 if s[0]=='F4'][0])<2]
prog_of={}; leader_hit={}
# programs are assigned by vertical order: labels top-to-bottom matched to F4 endpoints top-to-bottom (labels are black,
# so colour cannot carry identity); the leader-line endpoint is recorded separately to show which label it is nearest to
labs_sorted=labs.sort_values('cy').text.tolist(); ends_sorted=sorted(endF4,key=lambda e:e[1])
for (cc,_),lab in zip(ends_sorted,labs_sorted): prog_of[cc]=lab
for l in leaders:
    (xa,ya),(xb,yb)=sorted(l['pts'],key=lambda t:t[0]); cc=min(endF4,key=lambda e:abs(e[1]-ya))[0]
    d=labs.assign(dist=[max(0,t-yb,yb-bt) for t,bt in zip(labs.top,labs.bottom)]).sort_values('dist')
    leader_hit[prog_of[cc]]=f"nearest label {d.text.iloc[0]} ({d.dist.iloc[0]:.2f} pt), next {d.text.iloc[1]} ({d.dist.iloc[1]:.2f} pt)"
print("leader-line endpoints:",leader_hit)
print("colour->program via leader lines",prog_of)
for c in pts:
    cx=(c['x0']+c['x1'])/2; st=min(stx2,key=lambda s:abs(s[1]-cx))[0]; pr=prog_of.get(col(c['non_stroking_color']),'?'); v=np.polyval(b,cyc(c))
    if st=='F0': exp=0.0
    else:
        m=N[(N.program==pr)&(N.group==st)]; exp=m.est.iloc[0] if len(m) else np.nan
    out.append(dict(mark='nmf',key=f"{pr}:{st}",drawn=v,expected=exp))
O=pd.DataFrame(out); rng_b=abs(float(bot[0][1])-float(bot[-1][1]))
O['tol']=np.where(O.mark=='bar',np.maximum(0.5,0.005*rng_top),0.005*rng_b)
O['off']=(O.drawn-O.expected).abs()>O.tol
print(O.round(4).to_string()); print("marks checked",len(O),"off",int(O.off.sum()),"| printed top tick labels",signs_printed)
print("y-axis titles:",[ ''.join(W[(~W.upright)&(W.x0<16)&(W.top<60)].sort_values('top',ascending=False).text), ''.join(W[(~W.upright)&(W.x0<16)&(W.top>60)].sort_values('top',ascending=False).text)])
O.to_csv(f"out/fig3e_geometry_{tag}.tsv",sep="\t",index=False)
