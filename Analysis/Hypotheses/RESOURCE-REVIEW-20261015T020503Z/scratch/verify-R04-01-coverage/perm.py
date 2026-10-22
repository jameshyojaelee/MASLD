import sys; sys.path.insert(0, sys.argv[1]); from refit import *
N=int(sys.argv[2]); SEED=int(sys.argv[3])
rng=np.random.default_rng(SEED)
donors=el[el.donor.isin(d.donor.unique())].reset_index(drop=True)
# permuted stage matrix (N x donors), shuffled within dataset; column 0 = observed
P=np.tile(donors.stage.values,(N+1,1))
for ds,idx in donors.groupby("dataset").indices.items():
    for k in range(1,N+1): P[k,idx]=rng.permutation(donors.stage.values[idx])
pos={dn:i for i,dn in enumerate(donors.donor)}
T=np.zeros((117,N+1)); TH=np.zeros((117,N+1))
for j,(ct,mod) in enumerate(progs.itertuples(index=False)):
    sub=d[(d.cell_type==ct)&(d.module==mod)]
    ii=np.array([pos[x] for x in sub.donor]); y=sub.score.values; ds=sub.dataset.values
    S_=P[:,ii].copy(); yr=y.copy(); h0=np.zeros(len(y))
    for lev in set(ds):
        m=ds==lev; S_[:,m]-=S_[:,m].mean(1,keepdims=True); yr[m]-=yr[m].mean(); h0[m]=1/m.sum()
    n=len(y); rdf=n-(1+len(set(ds)))
    SS=(S_**2).sum(1); ok=SS>1e-10
    b=np.where(ok,(S_@yr)/np.where(ok,SS,1),0)
    e=yr[None,:]-b[:,None]*S_
    se=np.sqrt((e**2).sum(1)/rdf/np.where(ok,SS,1))
    t=np.where(ok,b/se,0)
    h=h0[None,:]+S_**2/np.where(ok,SS,1)[:,None]
    good=ok&np.all(1-h>1.5e-8,axis=1)
    v=((S_/np.where(ok,SS,1)[:,None])**2*(e/np.where(1-h>1.5e-8,1-h,1))**2).sum(1)
    th=np.where(good,b/np.sqrt(np.where(v>0,v,1)),0)
    T[j]=t; TH[j]=th
duc=np.where((progs.cell_type=="hepatocytes")&(progs.module==20))[0][0]
# check observed vs closed-form fit
sub=d[(d.cell_type=="hepatocytes")&(d.module==20)]; b0,p0,ph0,t0,th0=fit(sub)
print("observed t check:",T[duc,0],t0," hc3 t check:",TH[duc,0],th0)
for lab,M in (("OLS t",T),("HC3 t",TH)):
    obs=np.abs(M[:,0]); null=np.abs(M[:,1:])
    pp=(1+(null>=obs[:,None]).sum(1))/(N+1)
    q=bh(pp)
    mx=null.max(0)
    fwer=(1+(mx[None,:]>=obs[:,None]).sum(1))/(N+1)
    print(f"{lab}: Ductular perm p={pp[duc]:.5g} BH117 q={q[duc]:.5g} maxT FWER={fwer[duc]:.5g}; programs with BH117 q<0.05: {(q<0.05).sum()}; FWER<0.05: {(fwer<0.05).sum()}; floor 1/(N+1)={1/(N+1):.3g}")
    o=np.argsort(pp)[:4]
    for k in o: print("   ",progs.cell_type[k],progs.module[k],f"p={pp[k]:.4g} q={q[k]:.4g} fwer={fwer[k]:.4g}")
