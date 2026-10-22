import pandas as pd, numpy as np, sys
J="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3031/out"
S=sys.argv[1]
N={"MVP_ALT_EUR":406803,"MVP_ALT_AFR":115128,"MVP_ALT_AMR":56001,"MVP_ALT_EAS":6050,
   "MVP_AST_EUR":389192,"MVP_AST_AFR":112429,"MVP_AST_AMR":54250,"MVP_AST_EAS":5835,
   "BBJ_ALT":160000,"BBJ_AST":160000,"PanUKBB_AFR_ALT":6636,"PanUKBB_AFR_AST":6636,"PanUKBB_CSA_ALT":8876,"PanUKBB_CSA_AST":8876}
nk=pd.read_csv(f"{S}/null_keys.tsv",sep="\t",header=None,names=["chromosome","position"],dtype={"chromosome":str})
print("Implied phenotype SD = se*sqrt(2*af*(1-af)*N_registry), median over seeded null SNVs with 0.05<af<0.95")
for f,n in N.items():
    d=pd.read_csv(f"{J}/{f}.tsv",sep="\t",dtype={"chromosome":str})
    d=d.merge(nk,on=["chromosome","position"])
    af=pd.to_numeric(d.af,errors="coerce"); ok=(af>0.05)&(af<0.95)
    sdy=d.se[ok]*np.sqrt(2*af[ok]*(1-af[ok])*n)
    print(f"  {f:18s} n_snv={ok.sum():5d} median_sdy={np.median(sdy):.3f}")
# negative control: UKBB p>0.5 SNVs, oriented to UKBB allele1, slope through origin
def load(name):
    d=pd.read_csv(f"{J}/{name}.tsv",sep="\t",dtype={"chromosome":str})
    return d[(d.allele1.str.len()==1)&(d.allele2.str.len()==1)]
rng=np.random.default_rng(20261015)
print("Negative control: IVW slope of target beta on UKBB beta at UKBB p>0.5 SNVs (non-palindromic), 1,000 drawn (seed 20261015)")
for T in ["ALT","AST"]:
    u=load(f"UKBB_{T}").merge(nk,on=["chromosome","position"])
    u=u[(u.pval>0.5)&~u.apply(lambda r:{r.allele1,r.allele2} in ({"A","T"},{"C","G"}),axis=1)]
    for a in ["EUR","AFR","AMR","EAS"]:
        t=load(f"MVP_{T}_{a}")
        m=u.merge(t,on=["chromosome","position"],suffixes=("","_t"))
        same=(m.allele1==m.allele1_t)&(m.allele2==m.allele2_t); swap=(m.allele1==m.allele2_t)&(m.allele2==m.allele1_t)
        m=m[same|swap].copy(); m["by"]=np.where(swap[same|swap],-m.beta_t,m.beta_t)
        m=m.drop_duplicates(["chromosome","position"])
        m=m.iloc[rng.choice(len(m),size=min(1000,len(m)),replace=False)]
        w=1/m.se_t**2; x=m.beta.values; y=m.by.values
        b=np.sum(w*x*y)/np.sum(w*x*x); se=np.sqrt(1/np.sum(w*x*x))
        print(f"  {T} MVP_{a}: n={len(m)} slope={b:.4f} SE={se:.4f} P(|slope|<0.05 | true 0)={2*__import__('scipy').stats.norm.cdf(0.05/se)-1:.2f}")
