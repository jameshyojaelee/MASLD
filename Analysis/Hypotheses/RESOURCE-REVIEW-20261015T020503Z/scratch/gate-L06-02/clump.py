import pandas as pd, sys
S=sys.argv[1]
for T in ["ALT","AST"]:
    d=pd.read_csv(f"{S}/ukbb_{T}_leadrows.tsv",sep="\t",dtype={"chromosome":str})
    n0=len(d)
    snv=d[(d.allele1.str.len()==1)&(d.allele2.str.len()==1)].copy()
    gws=snv[snv.pval<5e-8].sort_values("pval")
    kept=[]
    for r in gws.itertuples():
        if any(k[0]==r.chromosome and abs(k[1]-r.position)<=500000 for k in kept): continue
        kept.append((r.chromosome,r.position))
    kd=gws.set_index(["chromosome","position"]).loc[kept].reset_index()
    pal=kd.apply(lambda r: {r.allele1,r.allele2} in ({"A","T"},{"C","G"}),axis=1)
    kd["palindromic"]=pal
    kd.to_csv(f"{S}/ukbb_{T}_clumped500kb.tsv",sep="\t",index=False)
    print(T,"rows",n0,"SNV",len(snv),"GWS SNV",len(gws),"clumped",len(kd),"palindromic",int(pal.sum()),"nonpal",int((~pal).sum()),"mean|beta|",round(kd.beta.abs().mean(),4),"median|beta|",round(kd.beta.abs().median(),4))
