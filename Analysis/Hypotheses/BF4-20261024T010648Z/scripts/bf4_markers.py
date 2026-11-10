#!/usr/bin/env python
"""BF4 step 1: ranked single-nucleus marker candidates per cell type, and the top-rate cell type of every gene
(PRESPEC_BF4.md). Final 30-gene sets are chosen in bf4_bulk.R after checking presence in all five cohorts."""
import sys

import numpy as np
import pandas as pd

OUT = sys.argv[1]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
BUILD = f"{ROOT}/Analysis/Hypotheses/BF3-20261023T154841Z/build-20261023T155101Z"
CG0 = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv"
DIAG = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"]
TYPES = ["Cholangiocytes", "Fibroblasts", "Endothelial cells", "Macrophages", "T cells", "Hepatocytes"]

z = np.load(f"{BUILD}/donor_type_counts.npz", allow_pickle=True)  # BF3 build output; genes saved as object array
counts, donors, types, genes = z["counts"], z["donors"].astype(str), z["types"].astype(str), z["genes"].astype(str)
dt = pd.read_csv(f"{BUILD}/donors.tsv", sep="\t").set_index("donor").loc[donors]
ok = (dt.n_nuclei >= 500) & dt.diagnosis.isin(DIAG) & dt.sex.isin(["Male", "Female"])
counts, dt = counts[ok.values], dt[ok]
nuc = pd.Series({t: dt[f"nuclei_{t}"].sum() for t in types})
big = [i for i, t in enumerate(types) if nuc[t] >= 1000]
btypes = list(types[big])
rate = counts[:, big, :].sum(0) / counts[:, big, :].sum((0, 2))[:, None]  # types x genes, per-UMI rate
print("donors", len(dt), "types", btypes, flush=True)

g = pd.read_csv(CG0, sep="\t")
g["signature_gene"] = g.signature_gene.fillna(False).astype(bool)
nmap = g.dropna(subset=["gene_name"]).groupby("gene_name").g.agg(list)
ens = {s: v[0] for s, v in nmap.items() if len(v) == 1}
sig = set(g.gene_name[g.signature_gene].dropna())
top10 = set(g.gene_name[g.cont_z.abs() >= g.cont_z.abs().quantile(0.9)].dropna())

R = pd.DataFrame(rate.T, index=genes, columns=btypes)
R = R[R.index.isin(ens.keys())]
rows = []
for t in TYPES:
    other = R.drop(columns=t).max(1)
    spec = R[t] / (other + 1e-7)
    c = pd.DataFrame({"cell_type": t, "gene": R.index, "ensembl": R.index.map(ens), "rate": R[t].values,
                      "specificity": spec.values, "signature_gene": R.index.isin(sig),
                      "top10pct_abs_cont_z": R.index.isin(top10)})
    c = c[(c.rate >= 1e-5) & ~c.signature_gene].sort_values("specificity", ascending=False).head(300)
    c["rank"] = np.arange(1, len(c) + 1)
    rows.append(c)
pd.concat(rows).to_csv(f"{OUT}/marker_candidates.tsv", sep="\t", index=False)

expressed = R[R.max(1) >= 1e-6]
src = pd.DataFrame({"gene": expressed.index, "ensembl": expressed.index.map(ens),
                    "top_type": expressed.idxmax(1).values,
                    "hep_ratio": (expressed["Hepatocytes"] / expressed.max(1)).values})
src.to_csv(f"{OUT}/gene_top_celltype.tsv", sep="\t", index=False)
print(src.top_type.value_counts(), flush=True)
