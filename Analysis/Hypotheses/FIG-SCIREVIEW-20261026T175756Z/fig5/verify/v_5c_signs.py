import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
C=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5C.tsv",sep="\t")
M=pd.read_csv(R+"figures/main/fig5_molecular_context/data/composite_mrna_protein_corrected.csv")[["gene","module"]]
C=C.merge(M,on="gene")
f=C[C.filled]
print(f.groupby(["module","feature"]).rho.apply(lambda s:f"{(s>0).sum()}+/{(s<0).sum()}-").unstack().to_string())
print("filled by module:",f.groupby("module").size().to_dict())
