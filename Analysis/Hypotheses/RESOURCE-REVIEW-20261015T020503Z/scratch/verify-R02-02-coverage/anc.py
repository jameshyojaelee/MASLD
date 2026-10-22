import pandas as pd
D="GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/"
m=pd.read_csv(D+"credible_set_members_annotated.tsv",sep="\t",dtype=str,keep_default_na=False)
print(m.shape)
g=m[m.pip_sum_gate_passed=="true"].copy()
print("gated rows",len(g))
g["unres"]=g.promoter_context_status=="unresolved_liftover"
print(g.groupby("ancestry").unres.agg(["sum","count","mean"]))
print(g.groupby("trait_scope").ancestry.value_counts())
s=g.groupby("study").unres.agg(["mean","count"]).sort_values("mean")
print(s.to_string())
print("consequence status:",g.consequence_annotation_status.value_counts().to_dict())
print(g.groupby("ancestry").consequence_annotation_status.value_counts())
print("abc by ancestry")
print(g.groupby("ancestry").abc_enhancer_overlap.value_counts().unstack(fill_value=0))
