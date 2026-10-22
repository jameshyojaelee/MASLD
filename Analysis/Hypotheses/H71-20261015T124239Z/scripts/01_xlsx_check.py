# H71 prespec: confirm the gate TSVs equal the two sheets of Kumak Supplementary Data 1.
import pandas as pd, sys
X="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/public_r4/biorxiv_2026.08.26.747338/media-1.xlsx"
G="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG10/"
xl=pd.ExcelFile(X); print("sheets",xl.sheet_names)
for n in xl.sheet_names:
    kx=xl.parse(n); print(n,kx.shape,list(kx.columns)[:4])
for s in ["NAS","Fibrosis"]:
    k=pd.read_csv(G+f"kumak_{s}.tsv",sep="\t")
    n=[x for x in xl.sheet_names if s.lower() in x.lower()]
    kx=xl.parse(n[0]); kx.columns=[str(c) for c in kx.columns]
    same=k.shape==kx.shape and (k.astype(str).values==kx.astype(str).values).all()
    print(s,"sheet",n[0],"tsv",k.shape,"xlsx",kx.shape,"cell-identical",same)
print(sys.version.split()[0],pd.__version__)
