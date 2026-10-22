"""Fig 3D geometry: recover PC1/PC2 of every drawn point from tick positions, colour class per point, and match the
drawn cloud to each recomputed PCA variant (sign flips allowed) with a linear assignment; report max displacement and
colour/stage agreement. Usage: python 11_fig3d_geometry.py <pdf> <tag>"""
import pdfplumber, numpy as np, pandas as pd, glob, sys, os
from scipy.optimize import linear_sum_assignment
pdf,tag=sys.argv[1],sys.argv[2]
O="out/"
p=pdfplumber.open(pdf).pages[0]
lines=p.lines; chars=p.chars
# axis ticks: short lines on the left (x ticks bottom) of main plot
yt=sorted([l for l in lines if abs(l['x1']-l['x0'])<2 and abs(l['top']-l['bottom'])<0.01 and l['x0']<30], key=lambda l:l['top'])
xt=sorted([l for l in lines if abs(l['top']-l['bottom'])<2 and abs(l['x1']-l['x0'])<0.01 and l['top']>130], key=lambda l:l['x0'])
def label_near(x=None,y=None):
    # tick labels: group chars by line
    cs=[c for c in chars if c['size']>5.9]
    if y is not None:
        row=[c for c in cs if abs((c['top']+c['bottom'])/2-y)<2.5 and c['x1']<27]
    else:
        row=[c for c in cs if abs((c['x0']+c['x1'])/2-x)<5 and 134<c['top']<142]
    return float(''.join(c['text'] for c in sorted(row,key=lambda c:c['x0'])))
ytv=[(l['top'],label_near(y=l['top'])) for l in yt]; xtv=[(l['x0'],label_near(x=l['x0'])) for l in xt]
print("y ticks",ytv); print("x ticks",xtv)
ay=np.polyfit([a for a,b in ytv],[b for a,b in ytv],1); ax=np.polyfit([a for a,b in xtv],[b for a,b in xtv],1)
print("residual tick fit y",np.max(np.abs(np.polyval(ay,[a for a,b in ytv])-[b for a,b in ytv])),"x",np.max(np.abs(np.polyval(ax,[a for a,b in xtv])-[b for a,b in xtv])))
pts=[c for c in p.curves if (c['x1']-c['x0'])<3 and c['stroke'] and c['top']>34 and c['x0']<127]
cols={(0.619608,0.619608,0.619608):'control',(0.788235,0.14902,0.368627):'F4',(0.835294,0.313725,0.478431):'F3',(0.882353,0.482353,0.588235):'F2',(0.929412,0.647059,0.698039):'F1',(0.980392,0.815686,0.807843):'F0'}
D=pd.DataFrame([dict(px=(c['x0']+c['x1'])/2,py=c['top'],col=cols.get(tuple(round(v,6) for v in c['stroking_color']),str(c['stroking_color']))) for c in pts])
D['PC1']=np.polyval(ax,D.px); D['PC2']=np.polyval(ay,D.py)
print("drawn points",len(D),D.col.value_counts().to_dict())
D.to_csv(O+f"fig3d_drawn_points_{tag}.tsv",sep="\t",index=False)
ros=pd.read_csv(O+"fig3d_recomputed_roster.tsv",sep="\t")
ros['expcol']=np.where(ros.arm=='strict_control','control','F'+ros.fibrosis_stage.fillna(-1).astype(int).astype(str))
print("expected colour classes:",ros.expcol.value_counts().to_dict())
best=None
rng=float(np.ptp(D.PC1))
for f in sorted(glob.glob(O+"pca_variants/"+os.environ.get("VARIANT_GLOB","*")+".tsv")):
    v=pd.read_csv(f,sep="\t").merge(ros,on='sample_id')
    for s1 in (1,-1):
        for s2 in (1,-1):
            A=np.c_[s1*v.PC1,s2*v.PC2]; B=D[['PC1','PC2']].values
            if len(A)!=len(B): continue
            C=np.sqrt(((A[:,None,:]-B[None,:,:])**2).sum(-1)); r,c=linear_sum_assignment(C)
            mx=C[r,c].max(); med=np.median(C[r,c])
            if best is None or med<best[1]: best=(os.path.basename(f),med,mx,s1,s2,r,c,v)
print("best variant",best[0],"median disp",round(best[1],3),"max disp",round(best[2],3),"signs",best[3],best[4])
_,med,mx,s1,s2,r,c,v=best
M=v.iloc[r].reset_index(drop=True).join(D.iloc[c].reset_index(drop=True),lsuffix='_x',rsuffix='_y')
M['disp']=np.sqrt((s1*M.PC1_x-M.PC1_y)**2+(s2*M.PC2_x-M.PC2_y)**2) if 'PC1_x' in M else np.nan
print("tolerance 0.5% of axis range =",round(0.005*rng,3))
M.to_csv(O+f"fig3d_matched_{tag}.tsv",sep="\t",index=False)
print("colour vs expected stage agreement:",(M.col==M.expcol).sum(),"of",len(M))
print(pd.crosstab(M.expcol,M.col))
print("points off (>0.5% range):",(M.disp>0.005*rng).sum())

# per-point within-colour check: points whose drawn colour disagrees with their own unit's stage
print("max disp",round(M.disp.max(),4),"median",round(M.disp.median(),4))
# recomputed variance
