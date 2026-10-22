import numpy as np, pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816"
d=pd.read_csv(R+"/multiplicity/corrected_locus_coordinates.tsv",sep="\t")
def cl(f):
    f=f.sort_values(["chrom","pos"]); ids=pd.Series(index=f.index,dtype=object)
    for c,b in f.groupby("chrom",sort=False):
        new=np.concatenate([[True],np.diff(b.pos.to_numpy())>1_000_000]); ids.loc[b.index]=[f"{c}:{i}" for i in np.cumsum(new)]
    return ids
d["pcl"]=cl(d)
e=d[d.ancestry=="EUR"].copy(); e["scl"]=cl(e)
ed=e[e.trait_scope=="tier1_direct_masld_pdff"]
print("EUR direct within-EUR clusters", ed.scl.nunique(), "-> pooled clusters", ed.pcl.nunique())
print(ed.groupby("scl").agg(pcl=("pcl","first"),pos=("pos","min"),n=("pos","size")).sort_values("pcl").to_string())
dd=d[d.trait_scope=="tier1_direct_masld_pdff"]
print("pooled direct clusters", dd.pcl.nunique(), "with EUR direct set", dd[dd.ancestry=="EUR"].pcl.nunique())
print("direct sets EUR share", (dd.ancestry=="EUR").mean())
