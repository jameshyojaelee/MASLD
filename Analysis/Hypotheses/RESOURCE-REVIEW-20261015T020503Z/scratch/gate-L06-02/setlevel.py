import pandas as pd, numpy as np, sys
S=sys.argv[1]; R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
J=f"{R}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3031/out"
m=pd.read_csv(f"{R}/figures/main/fig2_genetics/panels/Fig2E_ancestry_unique_coloc_GWS_membership.tsv",sep="\t")
c=pd.read_csv(f"{R}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gwas_name","gene","PP.H4.susie","top_snp","ancestry"])
e=c[(c.ancestry=="EUR")&c["PP.H4.susie"].notna()&c.gene.isin(m.gene)]
best=e.sort_values("PP.H4.susie",ascending=False).drop_duplicates("gene").merge(m[["gene","display_category"]],on="gene")
sp=best.top_snp.str.split(":",expand=True); best["chr"]=sp[0]; best["pos"]=pd.to_numeric(sp[1],errors="coerce")
def load(n):
    d=pd.read_csv(f"{J}/{n}.tsv",sep="\t",dtype={"chromosome":str}); return d[(d.allele1.str.len()==1)&(d.allele2.str.len()==1)]
def orient(ref,t):
    x=ref.merge(t,on=["chromosome","position"],suffixes=("","_t"))
    same=(x.allele1==x.allele1_t)&(x.allele2==x.allele2_t); swap=(x.allele1==x.allele2_t)&(x.allele2==x.allele1_t)
    x=x[same|swap].copy(); x["by"]=np.where(swap[same|swap],-x.beta_t,x.beta_t); x["sy"]=x.se_t
    return x.drop_duplicates(["chromosome","position"])
def ivw(x,y,s): w=1/s**2; return np.sum(w*x*y)/np.sum(w*x*x)
rng=np.random.default_rng(20261015)
for T in ["ALT","AST"]:
    cl=pd.read_csv(f"{S}/ukbb_{T}_clumped500kb.tsv",sep="\t",dtype={"chromosome":str}); cl=cl[~cl.palindromic]
    sub=best[best.gwas_name.isin({f"UKBB_{T}",f"MVP_{T}_EUR"})]
    tag={}
    for r in sub.itertuples():
        cc=cl[(cl.chromosome==r.chr)&((cl.position-r.pos).abs()<=500000)]
        for p in cc.position: tag.setdefault((r.chr,p),set()).add(r.display_category)
    lab=pd.DataFrame([(k[0],k[1],"EUR_only" if v=={"EUR_only"} else ("not_EUR_only" if "EUR_only" not in v else "mixed")) for k,v in tag.items()],columns=["chromosome","position","grp"])
    eur=orient(cl,load(f"MVP_{T}_EUR"))
    for a in ["AFR","AMR"]:
        an=orient(cl,load(f"MVP_{T}_{a}")).merge(eur[["chromosome","position","by","sy"]],on=["chromosome","position"],suffixes=("","_e")).merge(lab,on=["chromosome","position"])
        for g in ["EUR_only","not_EUR_only"]:
            k=an[an.grp==g]; n=len(k)
            if n<5: print(T,a,g,"n",n); continue
            x=k.beta.values; Rv=ivw(x,k.by.values,k.sy.values)/ivw(x,k.by_e.values,k.sy_e.values)
            bs=[]
            for i in range(5000):
                ix=rng.integers(0,n,n); bs.append(ivw(x[ix],k.by.values[ix],k.sy.values[ix])/ivw(x[ix],k.by_e.values[ix],k.sy_e.values[ix]))
            print(f"{T} MVP_{a} loci tagged by {g:13s} genes: n_loci={n:3d} R={Rv:.3f} boot_SE={np.std(bs,ddof=1):.3f}")
# PNPLA3/TM6SF2 per-variant ratios
for T in ["ALT","AST"]:
    e_=load(f"MVP_{T}_EUR")
    for a in ["AFR","AMR","EAS"]:
        t=load(f"MVP_{T}_{a}")
        for chrom,pos,nm in [("22",44324727,"PNPLA3 rs738409"),("19",19379549,"TM6SF2 rs58542926")]:
            x=e_[(e_.chromosome==chrom)&(e_.position==pos)].iloc[0]; y=t[(t.chromosome==chrom)&(t.position==pos)].iloc[0]
            r=y.beta/x.beta; se=abs(r)*np.sqrt((y.se/y.beta)**2+(x.se/x.beta)**2)
            print(f"{T} {nm} MVP_{a}/MVP_EUR per-SD ratio {r:.3f} (SE {se:.3f})")
