"""R11-01: (1) label-only participation ratio from the 99 GSE267145 grades; (2) small simulation of the
association-vector PR under one shared axis vs aspect-specific gene sets (5,000 genes, 3 reps, seed 20260923).
The full-size simulation is review job 4041 (30,605 genes); this is a light independent check."""
import numpy as np, pandas as pd
from scipy import stats
J="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"
OUT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F7b"
j=pd.read_csv(J,sep="\t"); L=j[["steatosis","ballooning","lobular_inflammation"]].astype(float).dropna()
print("participants",len(L))
def pr(M):
    ev=np.linalg.eigvalsh(M); return ev.sum()**2/(ev**2).sum(), ev[::-1]
R=stats.spearmanr(L).statistic; p,ev=pr(R)
print("label rho",R[0,1].round(4),R[0,2].round(4),R[1,2].round(4),"PR",round(p,4),"eig",ev.round(4))
print("C1 pass (all |rho|<0.80):",bool((np.abs(R[np.triu_indices(3,1)])<0.80).all()),"C2 pass (PR>=2):",p>=2.0)
rng=np.random.default_rng(20260923); n=len(L); G=5000
Lr=np.column_stack([stats.rankdata(L.iloc[:,k]) for k in range(3)]); Lz=(Lr-Lr.mean(0))/Lr.std(0)
def assoc_pr(X):
    Xr=np.apply_along_axis(stats.rankdata,0,X); Xr=(Xr-Xr.mean(0))/Xr.std(0)
    A=Xr.T@Lz/n          # genes x 3 Spearman associations
    return pr(stats.spearmanr(A).statistic)[0]
rows=[]
for model in ["noise","one_axis","separable"]:
  for frac in [0.05,0.2]:
    for beta in [0.5,1.0]:
      if model=="noise" and (frac,beta)!=(0.05,0.5): continue
      vals=[]
      for rep in range(3):
        X=rng.standard_normal((n,G)); k=int(frac*G)
        if model=="one_axis":
            s=Lz.sum(1); s=(s-s.mean())/s.std(); X[:,:k]+=beta*s[:,None]
        elif model=="separable":
            for a in range(3):
                idx=slice(a*k//3,(a+1)*k//3); X[:,idx]+=beta*Lz[:,[a]]
        vals.append(assoc_pr(X))
      rows.append(dict(model=model,frac=frac,beta=beta,pr_mean=np.mean(vals),pr_min=np.min(vals),pr_max=np.max(vals)))
t=pd.DataFrame(rows); print(t.round(4).to_string())
t.to_csv(f"{OUT}/r1101_separability_pr_small_simulation.tsv",sep="\t",index=False,float_format="%.4f")
pd.DataFrame([dict(statistic="label-only PR (99 participants)",value=p),dict(statistic="rho steatosis-ballooning",value=R[0,1]),dict(statistic="rho steatosis-lobular",value=R[0,2]),dict(statistic="rho ballooning-lobular",value=R[1,2])]).to_csv(f"{OUT}/r1101_label_pr.tsv",sep="\t",index=False,float_format="%.4f")
