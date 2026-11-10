import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/"
cs=pd.read_csv(R+"results/credible_sets.csv",usecols=["chromosome","position","locus","study"],dtype={"locus":str})
def nk(s): return s.astype(str).str.strip().str.rstrip("0").str.rstrip(".")
tot=[0,0,0]
for s in ["UKBB_ALT","UKBB_AST","UKBB_GGT"]:
    d=pd.read_csv(R+f"data/lead_snps/{s}_leadSNPs_tiered.tsv",sep="\t",dtype={"locus":str}); d["k"]=nk(d.locus)
    c=cs[cs.study==s].copy(); c["k"]=nk(c.locus); keys=set(c.k)
    nw=nk_=npos=0
    for ch,b in d.sort_values(["CHR","BP"]).groupby("CHR"):
        bp=b.BP.values; wid=np.cumsum(np.r_[True,np.diff(bp)>1e6])
        cp=c[c.chromosome.astype(str)==str(ch)].position.values
        for u in np.unique(wid):
            sel=wid==u; nw+=1
            nk_+=bool(set(b.k.values[sel])&keys)
            npos+=bool(((cp>=bp[sel].min()-5e5)&(cp<=bp[sel].max()+5e5)).any())
    print(s,nw,nk_,npos, "CS loci",c.locus.nunique()); tot=[tot[0]+nw,tot[1]+nk_,tot[2]+npos]
print(tot)
