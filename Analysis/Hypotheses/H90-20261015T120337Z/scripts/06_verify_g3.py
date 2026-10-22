"""H90 check: re-derive G3 coverage for Ductular injury, Stromal ECM and module 2 from the
raw GEO matrix files, using only the hepatocyte labels from job 6301. Reads no score."""
import gzip
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io

sys.path.insert(0, str(Path(__file__).parent))
import h90_common as C

obs = ad.read_h5ad(C.DERIVED / "gse289173_qc_nuclei.h5ad", backed="r").obs
hep = obs[obs.label.astype(str) == "Hepatocytes"]
hep_keys = set(hep.index)
feat = pd.read_csv(C.FEATURES_FL514, sep="\t", header=None, names=["id", "symbol", "type"])
mem = C.program_table()
rows = []
for ct, mod in (C.PRIMARY, C.SECONDARY, C.NC1):
    w = C.program_weights(mem, ct, mod)
    for g, wt in zip(w.mapped_symbol, w.original_l1_weight):
        rows.append(dict(module=mod, gene=g, weight=wt))
tab = pd.DataFrame(rows)
genes = sorted(set(tab.gene))
fidx = {g: np.flatnonzero(feat.symbol.to_numpy() == g) for g in genes}
det = {g: 0 for g in genes}
n_hep = 0
design = pd.read_csv(C.WRITE / "donor_design.tsv", sep="\t")
for r in design.itertuples():
    pre = f"{r.gsm}_{r.title}_"
    with gzip.open(C.RAW / f"{pre}barcodes.tsv.gz", "rt") as fh:
        bc = fh.read().split()
    keep = np.array([f"{r.gsm}_{b}" in hep_keys for b in bc])
    if not keep.any():
        continue
    m = scipy.io.mmread(C.RAW / f"{pre}matrix.mtx.gz").tocsr()[:, keep]
    n_hep += keep.sum()
    for g in genes:
        if len(fidx[g]):
            det[g] += int((np.asarray(m[fidx[g]].sum(axis=0)).ravel() > 0).sum())
tab["n_hep_nuclei"] = n_hep
tab["detected_nuclei"] = tab.gene.map(det)
tab["detection"] = tab.detected_nuclei / n_hep
tab["in_features"] = tab.gene.map(lambda g: len(fidx[g]) > 0)
tab["scored"] = tab.in_features & (tab.detection >= 0.01)
tab.to_csv(C.OUT / "G3_recheck_genes.tsv", sep="\t", index=False)
summ = tab.groupby("module").apply(lambda x: pd.Series({"n_scored": int(x.scored.sum()),
                                                         "retained_weight": x.weight[x.scored].sum(),
                                                         "n_mapped": len(x)}))
summ.to_csv(C.OUT / "G3_recheck_summary.tsv", sep="\t")
print(f"hepatocyte nuclei {n_hep}")
print(summ.to_string())
print(tab[tab.module == 20].sort_values("detection").to_string(index=False))
