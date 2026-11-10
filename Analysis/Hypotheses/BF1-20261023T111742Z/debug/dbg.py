import anndata as ad, numpy as np, scipy.sparse as sp
A = ad.read_h5ad("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/HRA_20.h5ad")
X = A.X
print(type(X), X.dtype, A.shape)
print(A.obs[["total_counts","n_genes"]].describe())
print("nan total", np.isnan(A.obs.total_counts).sum(), "nan ngenes", np.isnan(A.obs.n_genes.astype(float)).sum())
tot = np.asarray(X.sum(1)).ravel()
print("tot zero", (tot==0).sum(), "min", tot.min(), "dtype", tot.dtype)
print("tot vs obs total_counts max abs diff", np.abs(tot - A.obs.total_counts.values).max())
s = A.var.gene_symbol.astype(str).values
u, c = np.unique(s, return_counts=True); print("dup symbols", (c>1).sum())
