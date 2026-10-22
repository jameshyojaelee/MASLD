import pandas as pd, numpy as np, sys
S=sys.argv[1]; R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
v=pd.read_csv(f"{R}/GWAS/finemapping/results/credible_set_variant_consequences.csv")
v["chr"]=v.chr.astype(str)
print("annotation rows",len(v),"build: pos of PNPLA3 rs738409 in table:",v[(v.chr=="22")&(v.SYMBOL=="PNPLA3")][["pos","ref","alt","Consequence"]].head(3).to_string(index=False))
print(v["class"].value_counts().to_string())
for T in ["ALT","AST"]:
    cl=pd.read_csv(f"{S}/ukbb_{T}_clumped500kb.tsv",sep="\t",dtype={"chromosome":str})
    cl=cl[~cl.palindromic]
    m=cl.merge(v.drop_duplicates(["chr","pos"]),left_on=["chromosome","position"],right_on=["chr","pos"],how="left")
    print(T,"non-pal clumped leads",len(cl),"annotated",m["class"].notna().sum(),"class counts:",m["class"].value_counts().to_dict())
