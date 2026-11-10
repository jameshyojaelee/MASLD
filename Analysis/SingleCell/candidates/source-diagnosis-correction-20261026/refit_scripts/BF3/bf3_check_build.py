#!/usr/bin/env python
"""BF3 build check: recount one donor per dataset directly from atlas rows; print UMIs per nucleus and marker origin."""
import sys

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

BUILD = sys.argv[1]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
ATLAS = f"{ROOT}/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
PAIR = {ds: f"{ROOT}/data/{ds}/metadata/donor_pairing.csv" for ds in ["GSE202379", "GSE244832"]}

z = np.load(f"{BUILD}/donor_type_counts.npz", allow_pickle=True)  # own build output; genes saved as object array
counts, donors, types, genes = z["counts"], list(z["donors"]), list(z["types"]), z["genes"].astype(str)
dt = pd.read_csv(f"{BUILD}/donors.tsv", sep="\t").set_index("donor")

f = h5py.File(ATLAS, "r")
obs = read_elem(f["obs"])
R = f["raw"]["X"]
indptr = R["indptr"][:]
for ds in PAIR:
    d = dt[(dt.dataset == ds) & (dt.n_nuclei >= 500)].n_nuclei.idxmin()
    dp = pd.read_csv(PAIR[ds], dtype=str)
    runs = dp.loc[dp.donor_id == d.split("_", 1)[1], "rna_srrs"].iloc[0].split(";")
    rows = np.where(obs["sample"].astype(str).isin(runs).values)[0]
    direct = np.zeros((len(types), len(genes)))
    ct = obs.cell_type.astype(str).values
    for i in rows:
        lo, hi = indptr[i], indptr[i + 1]
        direct[types.index(ct[i]), R["indices"][lo:hi]] += R["data"][lo:hi]
    k = donors.index(d)
    print(ds, d, "rows", len(rows), "n_nuclei", dt.loc[d, "n_nuclei"],
          "max abs diff", np.abs(direct - counts[k]).max(), "total", direct.sum(), flush=True)
    assert len(rows) == dt.loc[d, "n_nuclei"] and np.array_equal(direct, counts[k])

umi = counts.sum((1, 2)) / dt.loc[donors, "n_nuclei"].values
print(pd.Series(umi, index=donors).groupby(dt.loc[donors, "dataset"].values).describe().round(0))
allc = counts.sum(0)
for gname in ["ALB", "APOA1", "KRT19", "COL1A1", "PTPRC", "CD163", "PECAM1"]:
    g = int(np.where(genes == gname)[0][0])
    sh = allc[:, g] / allc[:, g].sum()
    print(gname, types[int(sh.argmax())], round(float(sh.max()), 3))
