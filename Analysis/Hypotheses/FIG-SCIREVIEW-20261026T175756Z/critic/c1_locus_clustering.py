"""Critic check: are the 426 colocalized genes and the 26 overlap genes in Figure 1B
independent units, or do they share GWAS signals (same driving GWAS, top SNPs within 500 kb)?
Also: direct comparison of the colocalized list (26/426) with the expression-matched control (13/426)."""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
R = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
deg = pd.read_csv(R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
gl = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
deg["ens"] = deg.gene.str.split(".").str[0]
deg["DEG"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
U = deg.merge(gl, left_on="ens", right_on="ensembl")
U = U[U.susie_state_t12.isin(["positive", "tested_le_0.5"])].copy()
U["gen"] = U.susie_state_t12 == "positive"
print("U", len(U), "gen", U.gen.sum(), "DEG", U.DEG.sum(), "both", (U.gen & U.DEG).sum())

a = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
                usecols=["gwas_name", "gene", "ensembl", "PP.H4.susie", "top_snp"])
a["ensembl"] = a.ensembl.astype(str).str.split(".").str[0]
pos = U[U.gen][["ens", "susie_driving_gwas", "DEG"]]
m = pos.merge(a, left_on=["ens", "susie_driving_gwas"], right_on=["ensembl", "gwas_name"], how="left")
m = m.sort_values("PP.H4.susie", ascending=False).drop_duplicates("ens")
print("positives with driving row:", m.top_snp.notna().sum(), "of", len(m))
m["chr"] = m.top_snp.str.split(":").str[0]
m["bp"] = m.top_snp.str.split(":").str[1].astype(float)

# all positive pairs (any tier-1/2 GWAS with SuSiE PP.H4 > 0.5) to cluster by any shared signal
ap = a[(a["PP.H4.susie"] > 0.5) & a.ensembl.isin(pos.ens)].copy()
ap["chr"] = ap.top_snp.str.split(":").str[0]; ap["bp"] = ap.top_snp.str.split(":").str[1].astype(float)

def clusters(df, w=500_000):
    lab = {}; cid = 0
    for c, g in df.sort_values(["chr", "bp"]).groupby("chr"):
        last = -1e18
        for _, r in g.iterrows():
            if r.bp - last > w:
                cid += 1
            lab[r.ens] = cid; last = r.bp
    return lab

lab = clusters(m)
m["locus"] = m.ens.map(lab)
print("426 genes -> loci (driving top SNP chained within 500 kb):", m.locus.nunique())
sz = m.groupby("locus").size()
print("locus size distribution:", sz.value_counts().sort_index().to_dict())
ov = m[m.DEG]
print("26 overlap genes -> loci:", ov.locus.nunique())
osz = ov.groupby("locus").size()
print("overlap genes per locus:", osz[osz > 1].to_dict())
# genes in multi-gene loci that are DEG
print("overlap genes in loci with >1 colocalized gene:", int(ov.locus.map(sz).gt(1).sum()))
# per-locus DEG concordance: fraction of loci where >=1 DEG
locdeg = m.groupby("locus").DEG.agg(["sum", "size"])
print("loci with >=1 DEG:", int((locdeg["sum"] > 0).sum()), "of", len(locdeg))
# design effect: intra-locus correlation of DEG status among colocalized genes (ANOVA estimator)
k = locdeg["size"].values; y = locdeg["sum"].values; N = k.sum(); p = y.sum() / N
msb = ((y - k * p) ** 2 / k).sum() / (len(k) - 1)
msw = (y - y ** 2 / k).sum() / (N - len(k))
k0 = (N - (k ** 2).sum() / N) / (len(k) - 1)
icc = (msb - msw) / (msb + (k0 - 1) * msw)
deff = 1 + (((k ** 2).sum() / N) - 1) * icc
print(f"ICC of DEG status within locus {icc:.3f}; mean cluster-size term {((k**2).sum()/N):.2f}; design effect {deff:.2f}")

# matched control comparison
t = [[26, 426 - 26], [13, 426 - 13]]
print("Fisher colocalized 26/426 vs expression-matched control 13/426:", fisher_exact(t))
