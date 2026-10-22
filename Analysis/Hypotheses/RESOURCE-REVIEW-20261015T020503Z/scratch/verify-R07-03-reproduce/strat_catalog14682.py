import numpy as np, pandas as pd
CAT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet"
B = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
d = pd.read_parquet(CAT)
print(d.primary_evidence_class.value_counts().to_dict())
j = d[~d.primary_evidence_class.str.contains("untested", case=False)].copy()
j["g"] = j.ensembl_id.astype(str).str.split(".").str[0]
deg = pd.read_csv(B); deg["g"] = deg.gene.astype(str).str.split(".").str[0]; deg = deg.drop_duplicates("g")
j = j.merge(deg[["g", "AveExpr"]], on="g", how="left")
print("joint", len(j), "missing AveExpr", j.AveExpr.isna().sum())
cls = j.primary_evidence_class.str.lower()
j["gen"] = cls.str.contains("genetic|concordant"); j["de"] = cls.str.contains("state|concordant")
print("genetic", j.gen.sum(), "DE", j.de.sum(), "both", (j.gen & j.de).sum())
j = j.dropna(subset=["AveExpr"])
rng = np.random.default_rng(20261015)
for k in (5, 10, 20, 50):
    b = pd.qcut(j.AveExpr.rank(method="first"), k, labels=False)
    g = j.groupby(b).agg(n=("gen", "size"), G=("gen", "sum"), D=("de", "sum"))
    print(f"14,682 AveExpr {k} bins: expected {(g.G*g.D/g.n).sum():.2f}  (n used {len(j)})")
