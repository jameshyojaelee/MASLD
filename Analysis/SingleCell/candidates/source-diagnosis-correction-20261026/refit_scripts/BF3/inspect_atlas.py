import h5py, numpy as np
from anndata.io import read_elem
P="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
f=h5py.File(P,"r")
R=f["raw"]["X"]; print("raw X", dict(R.attrs), R["data"].shape, R["data"].dtype)
d=R["data"][:5000]; print("raw sample", d[:20], "frac int", np.mean(np.abs(d-np.round(d))<1e-6))
rv=read_elem(f["raw"]["var"]); print("raw var", rv.shape, rv.columns.tolist()[:10], rv.index[:5].tolist())
v=read_elem(f["var"]); print("var", v.shape, v.index[:5].tolist(), (v.index==rv.index).all() if len(v)==len(rv) else "len differ")
obs=read_elem(f["obs"]); print(obs.shape); print(obs.columns.tolist())
for c in obs.columns:
    if obs[c].dtype.name in("category","object","bool") and obs[c].nunique()<70: print(c, obs[c].value_counts().to_dict())
print(obs.groupby("dataset",observed=True)["sample"].nunique())
sub=obs[obs.dataset.isin(["GSE202379","GSE244832"])]
print("snRNA nuclei", len(sub)); 
for c in obs.columns:
    if "type" in c.lower() or "label" in c.lower() or "celltypist" in c.lower(): print("ct col",c, sub[c].nunique())
