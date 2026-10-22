"""V5 job 3 (a): donor pseudobulks by Ensembl ID, per-gene detection, and the H90
projection score (E3). Runs only after CHECKED_BEFORE_OUTCOME.txt says PROCEED: yes.
Writes into outputs/job3/. Follows PRESPEC_V5.md.

Pseudobulk sets (counts summed over cells, included donors only):
  hep            CellTypist hepatocytes (primary)
  hep_h90umi     CellTypist hepatocytes passing the H90 200-10,000 UMI rule (E7)
  hep_marker     marker-label hepatocytes (E8)
  <lineage>      each CellTypist lineage used by the 117 programs (E5, E6)
Usage: python -I 03_pseudobulk.py
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parent))
import v5_common as C

chk = C.RUN / "CHECKED_BEFORE_OUTCOME.txt"
if not chk.exists() or "PROCEED: yes" not in chk.read_text():
    sys.exit("STOP: CHECKED_BEFORE_OUTCOME.txt missing or not PROCEED: yes")
print(f"prespec sha256 {C.check_prespec()}", flush=True)
out = C.OUT / "job3"
out.mkdir(exist_ok=False)

A = ad.read_h5ad(C.DERIVED / "gse298719_qc_cells.h5ad")
AE = ad.read_h5ad(C.DERIVED / "gse298719_qc_cells_ensembl.h5ad")
if not (A.obs_names == AE.obs_names).all():
    sys.exit("STOP: symbol and Ensembl cell order differ")
obs = A.obs.reset_index(drop=True)
XS = A.X.tocsr()
XE = AE.X.tocsr()
genes = pd.Index(A.var_names.astype(str))
ens = pd.Index(AE.var_names.astype(str))
del A, AE
dq = pd.read_csv(C.OUT / "donor_qc.tsv", sep="\t", dtype={"id": str})
inc = sorted(dq.gsm[dq.included.astype(bool)])
inc_mask = obs.gsm.isin(inc).to_numpy()

sets = {
    "hep": (obs.label == "Hepatocytes").to_numpy(),
    "hep_h90umi": ((obs.label == "Hepatocytes") & (obs.umi <= 10000)).to_numpy(),
    "hep_marker": (obs.marker_lineage == "hepatocyte").to_numpy(),
}
for ct, lab in C.LINEAGE_LABEL.items():
    if ct != "hepatocytes":
        sets[ct] = (obs.label == lab).to_numpy()

ncell_rows, det_tabs = [], {}
for name, m in sets.items():
    m = m & inc_mask
    sub = obs[m]
    codes = pd.Categorical(sub.gsm, categories=inc).codes
    ind = sparse.csr_matrix((np.ones(m.sum()), (codes, np.arange(m.sum()))), shape=(len(inc), m.sum()))
    pb = (ind @ XE[m]).T.toarray()
    pd.DataFrame(pb, index=ens, columns=inc).rename_axis("gene_id").to_csv(out / f"pseudobulk_{name}.tsv.gz", sep="\t")
    det_tabs[name] = np.asarray((XE[m] > 0).mean(axis=0)).ravel() if m.sum() else np.zeros(len(ens))
    n = sub.groupby("gsm").size().reindex(inc).fillna(0).astype(int)
    ncell_rows += [dict(set=name, gsm=g, n_cells=int(v)) for g, v in n.items()]
pd.DataFrame(det_tabs, index=ens).rename_axis("gene_id").to_csv(out / "detection_by_set.tsv.gz", sep="\t")
pd.DataFrame(ncell_rows).to_csv(out / "cells_per_donor_set.tsv", sep="\t", index=False)

# ---- E3: H90 projection score of the three programs in CellTypist hepatocytes --------
m = sets["hep"] & inc_mask
hob = obs[m].reset_index(drop=True)
Xh = XS[m]
codes = pd.Categorical(hob.gsm, categories=inc).codes
gi = {g: i for i, g in enumerate(genes)}
det_h = np.asarray((Xh > 0).mean(axis=0)).ravel()
mem = C.membership()
rows = []
for uid in (C.PRIMARY, C.SECONDARY, C.NC1):
    w = C.weights(mem, uid)
    w = w[w.mapped_symbol.map(lambda g: g in gi and det_h[gi[g]] >= C.DET_MIN)]
    if len(w) < C.MIN_GENES:
        continue
    cols = [gi[g] for g in w.mapped_symbol]
    norm = C.log_normalize(Xh[:, cols].tocsr(), hob.umi.to_numpy(float)).astype(np.float64)
    n = norm.shape[0]
    mu = np.asarray(norm.mean(axis=0)).ravel()
    sq = norm.copy()
    sq.data **= 2
    sd = np.sqrt(np.clip((np.asarray(sq.mean(axis=0)).ravel() - mu ** 2) * n / (n - 1), 0, None))
    sd[sd <= 0] = 1.0
    ind = sparse.csr_matrix((np.ones(n), (codes, np.arange(n))), shape=(len(inc), n))
    dm = np.asarray((ind @ norm).todense()) / np.asarray(ind.sum(axis=1))
    z = (dm - mu) / sd
    wv = (w.original_l1_weight / w.original_l1_weight.sum()).to_numpy()
    rows += [dict(program_uid=uid, gsm=g, projection=s, n_genes=len(w)) for g, s in zip(inc, z @ wv)]
pd.DataFrame(rows).to_csv(out / "E3_projection_scores.tsv", sep="\t", index=False)
print(f"pseudobulks written for {len(inc)} donors: {', '.join(sets)}", flush=True)
