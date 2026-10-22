"""Fig 3F geometry: decode the indexed-colour rasters (heatmap 4x54, PP.H4 lane 1x54, two colourbars), invert each cell
colour to a value through its colourbar (piecewise-linear in RGB), read row labels and significance dots, and compare with
the recomputed matrix (stage log2FC clipped at +/-L, BH q<0.05 marks, best SuSiE PP.H4). Usage: <pdf> <tag>"""
import pdfplumber, numpy as np, pandas as pd, sys
pdf,tag=sys.argv[1],sys.argv[2]
p=pdfplumber.open(pdf).pages[0]
def decode(im):
    s=im['stream']; cs=s.attrs['ColorSpace']; pal=np.frombuffer(cs[3],dtype=np.uint8).reshape(-1,3)
    w,h=im['srcsize']; idx=np.frombuffer(s.get_data(),dtype=np.uint8)[:w*h].reshape(h,w)
    return pal[idx].astype(float)
ims=sorted(p.images,key=lambda i:(i['top'],i['x0']))
heat=[i for i in ims if i['srcsize'][1]>40 and i['srcsize'][0]>1][0]; lane=[i for i in ims if i['srcsize']==(1,heat['srcsize'][1])][0]
cbars=[i for i in ims if i['srcsize'][1]<40]; cb_h=max(cbars,key=lambda i:i['srcsize'][0]); cb_p=min(cbars,key=lambda i:i['srcsize'][0])
H=decode(heat); Lp=decode(lane)
# orientation: decide whether raster rows run top-down or bottom-up by testing both against the recomputed matrix below
FLIP=None; CH=decode(cb_h)[0]; CP=decode(cb_p)[0]
ticks=sorted([l for l in p.lines if l['top']>cb_h['bottom']-1],key=lambda l:l['x0'])
words=p.extract_words(); W=pd.DataFrame(words)
def tickval(x):
    c=W[(abs((W.x0+W.x1)/2-x)<6)&(W.top>ticks[0]['bottom'])]; return float(c.text.iloc[0].replace('−','-'))
tv=[(t['x0'],tickval(t['x0'])) for t in ticks]; print("colourbar ticks",tv)
def invert(rgb,cbar,img,t0,t1):
    # pixel centres of colourbar in page x; value linear in x between ticks t0 (x,val) and t1
    n=len(cbar); xs=img['x0']+(np.arange(n)+0.5)*(img['x1']-img['x0'])/n
    vals=t0[1]+(xs-t0[0])*(t1[1]-t0[1])/(t1[0]-t0[0])
    best=(1e9,None)
    for k in range(n-1):
        a,b=cbar[k],cbar[k+1]; d=b-a; tt=0 if (d@d)==0 else np.clip(((rgb-a)@d)/(d@d),0,1)
        e=np.linalg.norm(a+tt*d-rgb)
        if e<best[0]: best=(e,vals[k]+tt*(vals[k+1]-vals[k]))
    return best[1],best[0]
hv=[t for t in tv if t[0]<heat['x1']+5]; pv=[t for t in tv if t[0]>heat['x1']+5]
t0,t1=hv[0],hv[-1]; p0,p1=pv[0],pv[-1]
import os
if os.environ.get("HEAT_LIMIT"):
    L=float(os.environ["HEAT_LIMIT"]); t0=(t0[0],-L); t1=(t1[0],L); print("colourbar ends treated as +/-",L)
# row labels
labs=W[(W.x1<heat['x0'])&(W.x0>15)&(W.top<heat['bottom'])&(W.text.str.match(r'^[A-Z][A-Z0-9-]*$'))].copy(); labs['cy']=(labs.top+labs.bottom)/2
nr=heat['srcsize'][1]; rowh=(heat['bottom']-heat['top'])/nr
labs['row']=((labs.cy-heat['top'])/rowh).astype(int); labs=labs.sort_values('row')
genes=labs.text.tolist(); print("rows",len(genes))
ncol=heat['srcsize'][0]; colw=(heat['x1']-heat['x0'])/ncol
dots=[c for c in p.curves if (c['x1']-c['x0'])<3]
dotcells=set()
for c in dots:
    cx=(c['x0']+c['x1'])/2; cy=(min(q[1] for q in c['pts'])+max(q[1] for q in c['pts']))/2
    dotcells.add((int((cy-heat['top'])/rowh),int((cx-heat['x0'])/colw)))
R=pd.read_csv("out/fig3f_recomputed_matrix.tsv",sep="\t").set_index('gene_name')
contr=['F0_to_F1','F1_to_F2','F2_to_F3','F3_to_F4']
def build(Hm,Lm):
  rows=[]
  for i,g in enumerate(genes):
    for j,c in enumerate(contr):
        v,err=invert(Hm[i,j],CH,cb_h,t0,t1)
        exp=R.loc[g,c+'_lfc']; L=max(abs(t0[1]),abs(t1[1]))
        rows.append(dict(gene=g,contrast=c,drawn=v,expected_raw=exp,expected_clipped=float(np.clip(exp,-0.688,0.688)),fit_err=err,
            dot=(i,j) in dotcells,sig=R.loc[g,c+'_q']<0.05))
    v,err=invert(Lm[i,0],CP,cb_p,p0,p1)
    rows.append(dict(gene=g,contrast='PP.H4',drawn=v,expected_raw=R.loc[g,'pp4_best_susie'],expected_clipped=R.loc[g,'pp4_best_susie'],fit_err=err,dot=np.nan,sig=np.nan))
  return pd.DataFrame(rows)
cands={'top_down':build(H,Lp),'bottom_up':build(H[::-1],Lp[::-1])}
score={k:(v.drawn-v.expected_clipped).abs().median() for k,v in cands.items()}; print("orientation median abs diff",score)
FLIP=min(score,key=score.get); O=cands[FLIP]; print("raster orientation used:",FLIP)
tolH=0.5*(1/len(CH))*(t1[1]-t0[1])+0.005*(t1[1]-t0[1]); tolP=0.5/len(CP)+0.005
O['tol']=np.where(O.contrast=='PP.H4',tolP,tolH)
O['off']=(O.drawn-O.expected_clipped).abs()>O.tol
O['dot_mismatch']=(O.contrast!='PP.H4')&(O['dot'].astype(str)!=O['sig'].astype(str))
print(O[O.off|O.dot_mismatch].round(3).to_string())
print("cells checked",(O.contrast!='PP.H4').sum(),"off",O[(O.contrast!='PP.H4')].off.sum(),"| lane checked",(O.contrast=='PP.H4').sum(),"off",O[O.contrast=='PP.H4'].off.sum(),
      "| dots drawn",len(dots),"dot mismatches",O.dot_mismatch.sum(), "| genes recomputed but not drawn:",sorted(set(R.index)-set(genes)), "drawn not recomputed:",sorted(set(genes)-set(R.index)))
print("max |drawn-expected| heat",O[O.contrast!='PP.H4'].eval('drawn-expected_clipped').abs().max().round(4),"lane",O[O.contrast=='PP.H4'].eval('drawn-expected_clipped').abs().max().round(4))
O.to_csv(f"out/fig3f_geometry_{tag}.tsv",sep="\t",index=False)
pd.Series(genes).to_csv(f"out/fig3f_row_order_drawn_{tag}.tsv",sep="\t",index=False)
