"""audit-H58: fine-mapping window coverage using the lead file's own coordinates for windows and the
locus key (string, trailing zeros stripped on both sides) for the match, so no build mixing.
Locus strings are only compared, never parsed for a position. Counts only; no coding mass."""
import pandas as pd, numpy as np, os
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
reg=pd.read_csv(ROOT+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/snapshot/GWAS/finemapping/config/gwas_registry.tsv",sep="\t")
cs=pd.read_csv(ROOT+"GWAS/finemapping/results/credible_sets.csv",usecols=["study","locus"],dtype=str)
norm=lambda s: s.astype(str).str.strip().str.rstrip("0").str.rstrip(".")
cs["k"]=norm(cs.locus)
a=pd.read_csv(ROOT+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv",sep="\t",usecols=["study","trait_scope","ancestry"])
rows=[]
for _,r in a.drop_duplicates("study").iterrows():
    rr=reg[reg.study_name==r.study].iloc[0]; lp=ROOT+"GWAS/finemapping/"+str(rr.leadsnps_path)
    l=pd.read_csv(lp,sep="\t",dtype={"locus":str}) if os.path.exists(lp) else pd.DataFrame()
    keys=set(cs.loc[cs.study==r.study,"k"])
    if len(l)==0: rows.append((r.ancestry,r.trait_scope,r.study,0,0,len(keys),0)); continue
    l["k"]=norm(l.locus)
    nw=nf=0
    for ch,b in l.sort_values(["CHR","BP"]).groupby("CHR"):
        wid=np.cumsum(np.concatenate([[True],np.diff(b.BP.to_numpy())>1_000_000]))
        for w in np.unique(wid):
            nw+=1; nf+=int(bool(set(b.k.to_numpy()[wid==w]) & keys))
    rows.append((r.ancestry,r.trait_scope,r.study,nw,nf,len(keys),len(keys & set(l.k))))
d=pd.DataFrame(rows,columns=["ancestry","scope","study","lead_windows","windows_with_cs","cs_locus_keys","cs_keys_in_lead_file"])
pd.set_option("display.width",200)
print(d.sort_values(["ancestry","scope","study"]).to_string(index=False))
g=d.groupby(["ancestry","scope"])[["lead_windows","windows_with_cs"]].sum(); g["coverage"]=(g.windows_with_cs/g.lead_windows).round(3)
print(g.to_string())
