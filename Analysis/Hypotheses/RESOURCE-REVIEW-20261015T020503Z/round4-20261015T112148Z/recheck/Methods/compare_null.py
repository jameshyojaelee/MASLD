import pandas as pd, numpy as np, sys
O=sys.argv[1]
rep=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/Methods/repair/R11-04_R11-08_donor_relabel_null_repair.tsv.gz",sep="\t")
rows=[];ps={}
for L in ["hepatocyte","stellate"]:
  m=pd.read_csv(f"{O}/null_{L}.tsv",sep="\t")
  r=rep[rep.lineage==L].set_index("relabeling")
  m=m.set_index("j")
  agree={c:(m[a]==r[c].reindex(m.index)).mean() for a,c in [("same","same_sign_disc"),("sig244","n_sig_gse244832"),("supp","n_supported"),("disco","n_discordant")]}
  o=m[m.obs==1].iloc[0]
  n_same=(m.same>=o.same).sum(); n_supp=(m.supp>=o.supp).sum()
  ps[L]=n_same/len(m)
  rows.append(dict(lineage=L,n=len(m),obs_index=int(m.index[m.obs==1][0]),obs_same=int(o.same),obs_supp=int(o.supp),n_ge_same=int(n_same),p_same=n_same/len(m),
    null_mean_same=m.same.mean(),null_sd_same=m.same.std(),n_ge_supp=int(n_supp),p_supp=n_supp/len(m),null_mean_supp=m.supp.mean(),
    median_sig244=m.sig244.median(),**{f"agree_{k}":v for k,v in agree.items()}))
df=pd.DataFrame(rows)
p=np.array([ps["stellate"],ps["hepatocyte"],1.0]); o=np.argsort(p); q=np.empty(3); 
adj=np.minimum.accumulate((p[o]*3/np.arange(1,4))[::-1])[::-1]; q[o]=np.minimum(adj,1)
df["q_BH_three_rows"]=[q[1],q[0]]
p2=np.array([ps["hepatocyte"],ps["stellate"]]);o=np.argsort(p2);adj=np.minimum.accumulate((p2[o]*2/np.arange(1,3))[::-1])[::-1];q2=np.empty(2);q2[o]=adj
df["q_BH_two_rows"]=q2
df["binom_sd_indep"]=[np.sqrt(1016*.25),np.sqrt(15*.25)]
print(df.T.to_string()); df.to_csv(f"{O}/recheck_summary.tsv",sep="\t",index=False)
