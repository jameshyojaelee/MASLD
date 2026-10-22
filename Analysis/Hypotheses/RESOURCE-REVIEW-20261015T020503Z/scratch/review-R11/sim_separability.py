# Can the GSE267145 aspect-separability statistics distinguish separable from one-axis biology,
# given the observed label correlations? Uses the deposited 99-participant histology labels only.
import numpy as np, csv
from scipy.stats import rankdata
rng = np.random.default_rng(20261015)
J='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv'
rows=list(csv.DictReader(open(J),delimiter='\t'))
A=np.array([[float(r[k]) for k in ('steatosis','ballooning','lobular_inflammation')] for r in rows])
n=A.shape[0]; print('n participants',n)
R=np.column_stack([rankdata(A[:,k]) for k in range(3)])
Rz=(R-R.mean(0))/R.std(0)
lab=np.corrcoef(Rz.T); 
def pr(C):
    ev=np.linalg.eigvalsh(C); return ev.sum()**2/(ev**2).sum()
print('label Spearman matrix off-diag', np.round(lab[np.triu_indices(3,1)],4), 'label PR', round(pr(lab),4), 'closed form', round(9/(3+2*(lab[np.triu_indices(3,1)]**2).sum()),4))
G=20000
def assoc(X):
    # Spearman of each region (row) with each aspect
    Xr=np.apply_along_axis(rankdata,1,X); Xr=(Xr-Xr.mean(1,keepdims=True))/Xr.std(1,keepdims=True)
    return Xr@Rz/n  # G x 3
def stats(X):
    a=assoc(X)
    ar=np.column_stack([rankdata(a[:,k]) for k in range(3)])
    C=np.corrcoef(ar.T)
    return pr(C), C[np.triu_indices(3,1)]
z=(Rz)  # standardized rank labels as drivers
F=Rz.sum(1); F=(F-F.mean())/F.std()
out=[]
for scen in ['null','separable_disjoint','one_shared_axis']:
    for beta in ([0.0] if scen=='null' else [0.3,0.6,1.0]):
        prs=[];c3=[]
        for rep in range(3):
            E=rng.standard_normal((G,n))
            b=rng.standard_normal(G)*beta
            if scen=='null': X=E
            elif scen=='separable_disjoint':
                k=rng.integers(0,3,G); X=E+b[:,None]*z[:,k].T
            else: X=E+b[:,None]*F[None,:]
            p,c=stats(X); prs.append(p); c3.append(c)
        c3=np.mean(c3,0)
        print(f'{scen:20s} beta_sd={beta:.1f} PR mean={np.mean(prs):.4f} (reps {np.round(prs,4)}) assoc-vector pair Spearman s|b,s|i,b|i={np.round(c3,3)}')
