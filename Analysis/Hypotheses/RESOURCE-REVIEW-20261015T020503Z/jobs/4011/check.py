import numpy as np, pandas as pd
from scipy import stats
import pdfplumber
G="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/gse83452-intervention-response-20260831/results_20260831/participant_scores.tsv"
d=pd.read_csv(G,sep="\t")
print("rows",len(d), "arms", d.arm.value_counts().to_dict())
p2=d[(d.baseline_nash=="nash") & d.histology_delta.isin([-1.0,0.0])].copy()
p2["resolved"]=p2.histology_delta==-1.0
col="delta__gene_median_elastic_net"
print("P2 n",len(p2),"resolved",int(p2.resolved.sum()),"persisted",int((~p2.resolved).sum()))
print(pd.crosstab(p2.arm,p2.resolved))
z=stats.norm.ppf(0.975)+stats.norm.ppf(0.80)
def summ(x,lab):
    r=x[x.resolved][col].values; q=x[~x.resolved][col].values
    diff=r.mean()-q.mean(); sd=x[col].std(ddof=1)
    mde=2.802*sd*np.sqrt(1/len(r)+1/len(q)); mde_exact=z*sd*np.sqrt(1/len(r)+1/len(q))
    print(f"{lab}: n_res={len(r)} n_per={len(q)} diff={diff:.6f} sd={sd:.6f} MDE80(2.802)={mde:.6f} MDE80(exact z={z:.4f})={mde_exact:.6f}")
    return diff,mde,r,q
pd_,pm,_,_=summ(p2,"pooled")
dd,dm,r,q=summ(p2[p2.arm=="Diet"],"diet")
print(f"panel ratio diet diff / pooled MDE = {dd/pm:.4f}")
print(f"diet diff / diet MDE = {dd/dm:.4f}")
print(f"direction-oriented (correct sign = negative): -diff/diet MDE = {-dd/dm:.4f}")
rng=np.random.default_rng(1)
x=p2[p2.arm=="Diet"]; vals=x[col].values; res=x.resolved.values
B=10000; bs=[]
for i in range(B):
    ii=rng.integers(0,len(vals),len(vals))
    if res[ii].sum() and (~res[ii]).sum(): bs.append(vals[ii][res[ii]].mean()-vals[ii][~res[ii]].mean())
bs=np.array(bs)
print(f"own bootstrap (seed 1, {len(bs)} valid draws) 95% CI = [{np.percentile(bs,2.5):.4f}, {np.percentile(bs,97.5):.4f}]")
print(f"fraction of bootstrap draws <= -diet MDE ({-dm:.4f}): {(bs<=-dm).mean():.4f}; <= -pooled MDE: {(bs<=-pm).mean():.4f}; <= 0: {(bs<=0).mean():.4f}")
w=stats.ttest_ind(r,q,equal_var=False)
se=np.sqrt(r.var(ddof=1)/len(r)+q.var(ddof=1)/len(q))
dfw=se**4/((r.var(ddof=1)/len(r))**2/(len(r)-1)+(q.var(ddof=1)/len(q))**2/(len(q)-1))
t=stats.t.ppf(0.975,dfw)
print(f"Welch: t={w.statistic:.3f} p={w.pvalue:.4f} df={dfw:.2f} 95% CI=[{dd-t*se:.4f}, {dd+t*se:.4f}]")
# Rendered panel: point x positions
P="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig7_model_framework/fig7E_verdict_vocabulary_and_refusal_thresholds.pdf"
with pdfplumber.open(P) as pdf:
    pg=pdf.pages[0]
    words=pg.extract_words()
    ticks={w['text']:(w['x0']+w['x1'])/2 for w in words if w['text'] in ('0.0','0.5','1.0','1.5')}
    print("tick centers",ticks)
    a=(ticks['1.5']-ticks['0.0'])/1.5
    lab={}
    lines=pg.extract_text_lines()
    for L in lines: lab[L['text']]=(L['top']+L['bottom'])/2
    curves=[c for c in pg.curves]
    pts=[((c['x0']+c['x1'])/2,(c['top']+c['bottom'])/2) for c in curves if abs((c['x1']-c['x0'])-(c['bottom']-c['top']))<0.5 and (c['x1']-c['x0'])<6]
    print("n circle-like curves",len(pts))
    for k,yv in lab.items():
        if k.strip() in ("treatment response","signed eQTL direction","severity, female only","NAFLD endpoint"):
            near=min(pts,key=lambda p:abs(p[1]-yv)) if pts else None
            if near: print(f"{k.strip()}: point x -> ratio {(near[0]-ticks['0.0'])/a:.3f} (dy={near[1]-yv:.2f})")
