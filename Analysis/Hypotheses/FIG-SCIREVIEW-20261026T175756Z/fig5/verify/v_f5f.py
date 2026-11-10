import pandas as pd, numpy as np, itertools
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
F=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5F.tsv",sep="\t")
u=F[F.mark=="unit_point"]
for (prog,ds),g in u.groupby(["program","dataset"]):
    x=g.x_excess.values; obs=x.sum()
    if len(x)<=12:
        sums=[ (np.array(s)*np.abs(x)).sum() for s in itertools.product([1,-1],repeat=len(x))]
        p=np.mean(np.array(sums)>=obs-1e-15)
    print(prog,ds,"n",len(x),"pos",(x>0).sum(),"mean",round(x.mean(),4),"signflip one-sided P",p)
pp=F[F.mark=="pooled_point"][["program","dataset","x_excess","q_gene_null","null95_lo","null95_hi","signflip_p","printed_label"]]
print(pp.to_string())
