"""audit-H58: proposed gate G7 (fine-mapping completeness), threshold fixed at 0.80 before running.
Per study: GWS lead rows merged by single linkage (gap > 1 Mb starts a new window, same rule as
assign_clusters); a window counts as fine-mapped if credible_sets.csv has >= 1 row of that study on
that chromosome within [min lead BP - 500 kb, max lead BP + 500 kb]. Counts only; no coding mass."""
import pandas as pd, numpy as np, os
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
reg=pd.read_csv(ROOT+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/snapshot/GWAS/finemapping/config/gwas_registry.tsv",sep="\t")
cs=pd.read_csv(ROOT+"GWAS/finemapping/results/credible_sets.csv",usecols=["study","chromosome","position"])
a=pd.read_csv(ROOT+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv",sep="\t",usecols=["study","trait_scope","ancestry"])
st=a.drop_duplicates("study")
rows=[]
for _,r in st.iterrows():
    rr=reg[reg.study_name==r.study].iloc[0]; lp=ROOT+"GWAS/finemapping/"+str(rr.leadsnps_path)
    if not os.path.exists(lp): rows.append((r.ancestry,r.trait_scope,r.study,np.nan,np.nan)); continue
    l=pd.read_csv(lp,sep="\t")
    if len(l)==0: rows.append((r.ancestry,r.trait_scope,r.study,0,0)); continue
    c=cs[cs.study==r.study]
    nw=nf=0
    for ch,b in l.sort_values(["CHR","BP"]).groupby("CHR"):
        bp=b.BP.to_numpy(); new=np.concatenate([[True],np.diff(bp)>1_000_000]); wid=np.cumsum(new)
        cc=c[c.chromosome.astype(str)==str(ch)].position.to_numpy()
        for w in np.unique(wid):
            lo,hi=bp[wid==w].min()-500_000,bp[wid==w].max()+500_000
            nw+=1; nf+=int(((cc>=lo)&(cc<=hi)).any())
    rows.append((r.ancestry,r.trait_scope,r.study,nw,nf))
d=pd.DataFrame(rows,columns=["ancestry","scope","study","windows","finemapped"])
pd.set_option("display.width",200)
print(d.sort_values(["ancestry","scope","study"]).to_string(index=False))
g=d.groupby(["ancestry","scope"])[["windows","finemapped"]].sum(); g["coverage"]=g.finemapped/g.windows
print(g.to_string())
