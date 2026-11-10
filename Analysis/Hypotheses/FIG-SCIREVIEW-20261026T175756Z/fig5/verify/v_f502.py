import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
d=pd.read_csv(R+"figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_pair_rows_with_background.tsv",sep="\t")
g=d[["gene","ensembl"]].drop_duplicates()
print("D genes (ensembl):",g.ensembl.nunique(),"symbols:",g.gene.nunique())
U=set(open(R+"Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt").read().split())
print("universe size",len(U))
B=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5B.tsv",sep="\t")
Bp=B[B.mark=="scatter_point"] if "mark" in B else B
print("5B rows",len(B),"marks",B.mark.value_counts().to_dict())
ins=g.gene.isin(U); print("D symbols in universe:",ins.sum(),"not:",(~ins).sum())
print("D in 5B by ensembl:",g.ensembl.isin(set(B.ensembl_id)).sum(),"by symbol:",g.gene.isin(set(B.gene)).sum())
print("not-in-universe examples:",sorted(g.gene[~ins])[:20])
# genetic universe
UG=set(open(R+"Analysis/Spatial/results/universe_validation/universe_genetic.txt").read().split())
print("universe_genetic size",len(UG),"D in genetic universe",g.gene.isin(UG).sum())
E=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5E.tsv",sep="\t")
print("5E rows",len(E),"display",E.display.sum(), E[E.display].block.value_counts().to_dict())
print(E[E.display].groupby("block").robust.value_counts().to_dict())
print("display_rule",E.display_rule.unique())
# Is any lower-block program open (n_open_genes>0)?
print("lower block n_open>0:",(E[(E.display)&(E.block=="protein_robust_only")].n_open_genes>0).sum())
# upper block: are all with n_open_genes>0 displayed?
print("programs with n_open_genes>0:",(E.n_open_genes>0).sum(),"displayed:",E[E.n_open_genes>0].display.sum())
print("robust total",E.robust.sum(),"robust displayed",E[E.robust==True].display.sum())
