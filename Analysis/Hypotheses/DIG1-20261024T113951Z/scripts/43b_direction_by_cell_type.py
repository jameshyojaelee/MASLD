"""DIG1 job 43b (exploratory, post hoc; follows 41c and 43a). Does the genetic-vs-disease direction
test behave differently for genes made mainly by hepatocytes? Motivation: in §29 several discordant
genes (GAS6, F2RL1) rise in disease and may be made by non-parenchymal cells [Speculation].
Genetic sign: 41c gene table (SuSiE top pair, tier 1-2, one gene per 1 Mb locus).
Hepatocyte enrichment: GSE298719 Flex, 13 donors. Per gene, donor-mean log1p CP10k per cell type
(cell types with >= 20 cells in a donor), averaged over donors; hep_ratio = hepatocyte mean minus the
highest other cell type mean. Hepatocyte-enriched = hep_ratio > 0 (fixed before outcomes).
Tests (descriptive, Fisher exact; no Holm family since the split is post hoc):
  S1  concordance OR among hepatocyte-enriched genes (disease sign from padj < 0.05, unadjusted fit).
  S2  concordance OR among other genes.
  S3  2x2 concordant x enriched (are hepatocyte-enriched genes more often concordant?).
  S4  as S1-S3 with the composition-adjusted disease sign (43a, bh_adj < 0.05).
Usage: python -I 43b_direction_by_cell_type.py <new_out_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)

G = pd.read_csv(RUN / "out_direction_posthoc_controls/gene_directions_susie_top_dedup.tsv", sep="\t")
C = pd.read_csv(RUN / "out_bulk_composition_adjusted/gene_disease_composition_adjusted.tsv.gz", sep="\t")
G = G.merge(C[["ensembl", "beta_unadj", "bh_unadj", "beta_adj", "bh_adj"]], on="ensembl", how="left")

A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included]
A = A[A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
genes = pd.Index(A.var_names)
want = [x for x in G.gene.astype(str) if x in genes]
X = sparse.csr_matrix(A.X).astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
Ls = L[:, [genes.get_loc(x) for x in want]].toarray()
lab = A.obs.label.astype(str).to_numpy(); gsm = A.obs.gsm.astype(str).to_numpy()
rows = []
for d in dq.gsm:
    for t in np.unique(lab[gsm == d]):
        mt = (gsm == d) & (lab == t)
        if mt.sum() >= 20:
            rows.append(pd.Series(Ls[mt].mean(0), index=want, name=(d, t)))
M = pd.DataFrame(rows); M.index = pd.MultiIndex.from_tuples(M.index, names=["gsm", "label"])
ctm = M.groupby(level="label").mean().T
ctm.to_csv(OUT / "celltype_mean_by_gene.tsv", sep="\t")
oth = ctm.drop(columns="Hepatocytes")
spec = pd.DataFrame(dict(gene=ctm.index, hep_mean=ctm.Hepatocytes, max_other=oth.max(1), top_other=oth.idxmax(1)))
spec["hep_ratio"] = spec.hep_mean - spec.max_other
spec["hep_enriched"] = spec.hep_ratio > 0
G = G.merge(spec, on="gene", how="left")
G.to_csv(OUT / "gene_directions_with_celltype.tsv", sep="\t", index=False)


def fisher(x, dcol):
    ds = np.sign(x[dcol]); gs = x.gen_sign
    a = int(((gs > 0) & (ds > 0)).sum()); b = int(((gs > 0) & (ds < 0)).sum()); c = int(((gs < 0) & (ds > 0)).sum()); d = int(((gs < 0) & (ds < 0)).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]]) if len(x) else (np.nan, np.nan)
    return dict(n=len(x), frac_concordant=(a + d) / len(x) if len(x) else np.nan, odds_ratio=orr, p_fisher=p)


res = []
for lab_, dcol, filt in [("unadjusted padj<0.05", "logFC", lambda x: x.padj < 0.05),
                         ("composition-adjusted bh_adj<0.05", "beta_adj", lambda x: x.bh_adj < 0.05)]:
    x = G[filt(G) & G.hep_enriched.notna()].copy()
    x["conc"] = np.sign(x[dcol]) * x.gen_sign > 0
    for nm, sub in [("hepatocyte-enriched", x[x.hep_enriched.astype(bool)]), ("other", x[~x.hep_enriched.astype(bool)]), ("all", x)]:
        res.append(dict(disease=lab_, stratum=nm, **fisher(sub, dcol)))
    t = pd.crosstab(x.hep_enriched.astype(bool), x.conc).reindex(index=[True, False], columns=[True, False], fill_value=0)
    orr, p = stats.fisher_exact(t.to_numpy())
    res.append(dict(disease=lab_, stratum="S3 enriched x concordant", n=int(t.to_numpy().sum()), frac_concordant=np.nan, odds_ratio=orr, p_fisher=p))
R = pd.DataFrame(res); R.to_csv(OUT / "concordance_by_celltype.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 300)
print(R.round(4).to_string())
x = G[G.padj < 0.05].copy(); x["conc"] = np.sign(x.logFC) * x.gen_sign > 0
print(x.sort_values("hep_ratio")[["gene", "gwas", "gen_sign", "logFC", "beta_adj", "bh_adj", "hep_mean", "max_other", "top_other", "hep_ratio", "conc"]].round(3).to_string())
