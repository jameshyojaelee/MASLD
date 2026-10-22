import pandas as pd, numpy as np, sys
from scipy import stats
J="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3031/out"
S=sys.argv[1]
N={"MVP_ALT_EUR":406803,"MVP_ALT_AFR":115128,"MVP_ALT_AMR":56001,"MVP_ALT_EAS":6050,
   "MVP_AST_EUR":389192,"MVP_AST_AFR":112429,"MVP_AST_AMR":54250,"MVP_AST_EAS":5835,
   "BBJ_ALT":160000,"BBJ_AST":160000,"PanUKBB_AFR_ALT":6636,"PanUKBB_AFR_AST":6636,"PanUKBB_CSA_ALT":8876,"PanUKBB_CSA_AST":8876}
comp={"A":"T","T":"A","C":"G","G":"C"}
def load(name):
    d=pd.read_csv(f"{J}/{name}.tsv",sep="\t",dtype={"chromosome":str,"allele1":str,"allele2":str})
    d=d[(d.allele1.str.len()==1)&(d.allele2.str.len()==1)]
    return d
def orient(ref, tgt):
    # ref: UKBB rows (allele1 = effect). tgt: other file. match on chr:pos + unordered alleles, no strand flips.
    m=ref.merge(tgt,on=["chromosome","position"],suffixes=("","_t"))
    same=(m.allele1==m.allele1_t)&(m.allele2==m.allele2_t)
    swap=(m.allele1==m.allele2_t)&(m.allele2==m.allele1_t)
    m=m[same|swap].copy()
    m["by"]=np.where(swap[same|swap], -m.beta_t, m.beta_t)
    m["sy"]=m.se_t
    if "af" in m.columns:
        af=pd.to_numeric(m.af,errors="coerce")
        m["eaf_ukbbA1"]=np.where(swap[same|swap],1-af,af)
    return m.drop_duplicates(["chromosome","position"])
def ivw(x,y,sy):
    w=1/sy**2; sxx=np.sum(w*x*x); b=np.sum(w*x*y)/sxx; se=np.sqrt(1/sxx)
    phi=np.sum(w*(y-b*x)**2)/(len(x)-1)
    return b,se,phi
rng=np.random.default_rng(20261015)
out=[]
for T in ["ALT","AST"]:
    cl=pd.read_csv(f"{S}/ukbb_{T}_clumped500kb.tsv",sep="\t",dtype={"chromosome":str})
    cl=cl[~cl.palindromic]
    eur=orient(cl,load(f"MVP_{T}_EUR"))
    res={}
    for tgt in [f"MVP_{T}_EUR",f"MVP_{T}_AFR",f"MVP_{T}_AMR",f"MVP_{T}_EAS",f"BBJ_{T}",f"PanUKBB_AFR_{T}",f"PanUKBB_CSA_{T}"]:
        m=orient(cl,load(tgt))
        b,se,phi=ivw(m.beta.values,m.by.values,m.sy.values)
        res[tgt]=m
        # paired ratio vs MVP EUR on shared loci, locus bootstrap
        k=m.merge(eur[["chromosome","position","by","sy"]],on=["chromosome","position"],suffixes=("","_eur"))
        x=k.beta.values; ya=k.by.values; sa=k.sy.values; ye=k.by_eur.values; se_e=k.sy_eur.values
        ba,_,_=ivw(x,ya,sa); be,_,_=ivw(x,ye,se_e); R=ba/be
        B=10000; n=len(k); Rs=np.empty(B)
        for i in range(B):
            idx=rng.integers(0,n,n)
            Rs[i]=ivw(x[idx],ya[idx],sa[idx])[0]/ivw(x[idx],ye[idx],se_e[idx])[0]
        seR=Rs.std(ddof=1)
        zH=stats.norm.ppf(1-0.05/6/2)+stats.norm.ppf(0.8); z05=stats.norm.ppf(0.975)+stats.norm.ppf(0.8)
        out.append(dict(trait=T,target=tgt,n_loci=len(m),slope_on_ukbb=round(b,4),se_model=round(se,4),overdisp_phi=round(phi,2),
                        se_x_sqrtphi=round(se*np.sqrt(max(phi,1)),4),n_paired=n,R_vs_MVP_EUR=round(R,4),R_boot_se=round(seR,4),
                        R_ci95=f"{np.percentile(Rs,2.5):.3f}-{np.percentile(Rs,97.5):.3f}",
                        MDE_R_alpha05=round(z05*seR,3),MDE_R_holm6=round(zH*seR,3)))
    # per-locus deviations for MVP ancestries
    for a in ["AFR","AMR","EAS"]:
        m=res[f"MVP_{T}_{a}"].merge(eur[["chromosome","position","by","sy"]],on=["chromosome","position"],suffixes=("","_eur"))
        Rh=[o for o in out if o["target"]==f"MVP_{T}_{a}"][0]["R_vs_MVP_EUR"]
        d=m.by-Rh*m.by_eur; sd=np.sqrt(m.sy**2+Rh**2*m.sy_eur**2)
        p=2*stats.norm.sf(np.abs(d/sd)); q=stats.false_discovery_control(p,method="bh")
        r=m.by/m.by_eur; ser=np.abs(r)*np.sqrt((m.sy/m.by)**2+(m.sy_eur/m.by_eur)**2)
        lo=r-1.645*ser; hi=r+1.645*ser
        equiv=((lo>0.5)&(hi<1.5)).sum()
        print(f"{T} MVP_{a}: loci {len(m)}; per-locus deviation BH q<0.05: {(q<0.05).sum()}; median SE of deviation {np.median(sd):.4f}; median per-locus ratio SE {np.median(ser):.3f}; loci with 90% CI of ratio inside 0.5-1.5 (TOST equivalence): {equiv}")
pd.DataFrame(out).to_csv(f"{S}/transfer_ratios.tsv",sep="\t",index=False)
pd.set_option("display.width",250); pd.set_option("display.max_columns",30)
print(pd.DataFrame(out).to_string(index=False))
