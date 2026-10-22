import pandas as pd, numpy as np
from pathlib import Path
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
D=B/"Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
srr2donor={}
for ds in ("GSE244832","GSE185477","GSE202379","GSE136103"):
    dp=pd.read_csv(B/f"data/{ds}/metadata/donor_pairing.csv",dtype=str)
    print(ds, dp.shape, list(dp.columns)[:8])
    for _,r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s!="nan": srr2donor[s.strip()]=f"{ds}_{r['donor_id']}"
meta=pd.read_csv(D/"sample_celltype_metadata.tsv",sep="\t")
print(meta.dataset.value_counts())
for ln in ["Hepatocytes","Cholangiocytes","Fibroblasts","Endothelial_cells","Macrophages","T_cells"]:
    cols=pd.read_csv(D/f"pseudobulk_{ln}.tsv.gz",sep="\t",index_col=0,nrows=1).columns
    mapped=[srr2donor.get(c) for c in cols]
    m=meta[meta.cell_type==ln.replace("_"," ")].set_index("sample")
    ds=[m.dataset.get(c, "NA") if c in m.index else "NA" for c in cols]
    df=pd.DataFrame({"col":cols,"donor":[d if d else c for d,c in zip(mapped,cols)],"mapped":[d is not None for d in mapped],"ds":ds})
    print(ln, "cols",len(cols),"donors",df.donor.nunique(),"unmapped cols",(~df.mapped).sum())
    print("  donors by ds:", df.groupby("ds").donor.nunique().to_dict(), " cells min/median:", (m.n_cells.min() if len(m) else None), (m.n_cells.median() if len(m) else None))
