"""Fig 4B recomputation from the integrated atlas: cell-type proportions (16 labels), per-type UMAP centroids, and the
cells/donors per coarse stage that a Steatohepatitis/Healthy donor-averaged density would use (2026-10-14 corrected roster)."""
import h5py, numpy as np, pandas as pd
root = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out = f"{root}/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
f = h5py.File("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad", "r")
def cat(k):
    g = f["obs"][k]
    if isinstance(g, h5py.Group):
        return pd.Categorical.from_codes(g["codes"][:], [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
    return pd.Series(g[:]).str.decode("utf-8") if g.dtype.kind == "S" or g.dtype.kind == "O" else g[:]
ct = pd.Series(cat("cell_type")); ds = pd.Series(cat("dataset")); smp = pd.Series(cat("sample"))
U = f["obsm"]["X_umap"][:]
T = pd.DataFrame(dict(cell_type=ct.astype(str), dataset=ds.astype(str), sample=smp.astype(str), u1=U[:, 0], u2=U[:, 1]))
prop = T.cell_type.value_counts().rename("n_cells").to_frame(); prop["frac"] = prop.n_cells / len(T)
cen = T.groupby("cell_type")[["u1", "u2"]].median()
prop = prop.join(cen); print(prop.to_string())
prop.to_csv(f"{out}/r4b_celltype_props.tsv", sep="\t")
M = pd.read_csv(f"{root}/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv", sep="\t")
print(M.columns.tolist()[:30])
print(M.disease_stage_coarse.value_counts().to_dict())
