#!/usr/bin/env python
"""BF3 step 1: donor x cell-type pseudobulk counts for GSE202379 and GSE244832 (PRESPEC_BF3.md)."""
import hashlib
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

OUT = sys.argv[1]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
ATLAS = f"{ROOT}/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
PAIR = {ds: f"{ROOT}/data/{ds}/metadata/donor_pairing.csv" for ds in ["GSE202379", "GSE244832"]}
DIAG = (f"{ROOT}/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/"
        "donor_metadata_stage_corrected_donor.tsv")
SEX = f"{ROOT}/Analysis/SingleCell/metadata/snrna_donor_sex.csv"
CHUNK = 100_000


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


with open(f"{OUT}/input_sha256.txt", "w") as f:
    for p in [*PAIR.values(), DIAG, SEX]:
        f.write(f"{sha(p)}  {p}\n")
    f.write(f"not hashed (40 GB)  {ATLAS}\n")

f = h5py.File(ATLAS, "r")
obs = read_elem(f["obs"])
genes = np.asarray(read_elem(f["raw"]["var"]).index.astype(str))
keep = obs.dataset.isin(list(PAIR)).values
sel = obs[keep]
groups = pd.MultiIndex.from_frame(sel[["sample", "cell_type"]].astype(str)).unique()
gidx = pd.Series(np.arange(len(groups)), index=groups)
row_group = np.full(len(obs), -1, dtype=np.int64)
row_group[keep] = groups.get_indexer(pd.MultiIndex.from_frame(sel[["sample", "cell_type"]].astype(str)))
assert (row_group[keep] >= 0).all()
print("nuclei", keep.sum(), "run x type groups", len(groups), flush=True)

R = f["raw"]["X"]
indptr = R["indptr"][:]
G = len(genes)
acc = np.zeros((len(groups), G))
for a in range(0, len(obs), CHUNK):
    b = min(a + CHUNK, len(obs))
    rg = row_group[a:b]
    if (rg >= 0).sum() == 0:
        continue
    lo, hi = indptr[a], indptr[b]
    X = sp.csr_matrix((R["data"][lo:hi], R["indices"][lo:hi], indptr[a:b + 1] - lo), shape=(b - a, G))
    m = rg >= 0
    M = sp.csr_matrix((np.ones(m.sum()), (rg[m], np.where(m)[0])), shape=(len(groups), b - a))
    acc += (M @ X).toarray()
    print(a, flush=True)
f.close()
assert np.allclose(acc, np.round(acc))

# run -> donor
rmap = {}
for ds, p in PAIR.items():
    dp = pd.read_csv(p, dtype=str)
    for d, srrs in zip(dp.donor_id, dp.rna_srrs):
        for s in srrs.split(";"):
            assert s not in rmap, s
            rmap[s] = f"{ds}_{d}"
runs = sorted(sel["sample"].astype(str).unique())
missing = [r for r in runs if r not in rmap]
assert not missing, missing
types = sorted(sel.cell_type.astype(str).unique())
donors = sorted({rmap[r] for r in runs})
di, ti = {d: i for i, d in enumerate(donors)}, {t: i for i, t in enumerate(types)}
counts = np.zeros((len(donors), len(types), G), dtype=np.float64)
for (run, ct), k in gidx.items():
    counts[di[rmap[run]], ti[ct]] += acc[k]
nuc = sel.assign(donor=sel["sample"].astype(str).map(rmap)).groupby(["donor", "cell_type"], observed=True).size()
nuc = nuc.unstack(fill_value=0).reindex(index=donors, columns=types, fill_value=0)

diag = pd.read_csv(DIAG, sep="\t").set_index("sample")
sex = pd.read_csv(SEX)
sex = sex[sex["sample"].isin(runs)].assign(donor=lambda d: d["sample"].map(rmap))
sexd = sex.groupby("donor").inferred_sex.agg(lambda s: s.mode().iloc[0] if s.value_counts().iloc[0] > len(s) / 2 else "tie")
dt = pd.DataFrame({"donor": donors, "dataset": [d.split("_")[0] for d in donors],
                   "n_runs": [sum(rmap[r] == d for r in runs) for d in donors],
                   "n_nuclei": nuc.sum(1).values,
                   "diagnosis": diag.reindex(donors).disease_stage_coarse.values,
                   "sex": sexd.reindex(donors).values})
dt = pd.concat([dt, nuc.add_prefix("nuclei_").reset_index(drop=True)], axis=1)
dt.to_csv(f"{OUT}/donors.tsv", sep="\t", index=False)
np.savez_compressed(f"{OUT}/donor_type_counts.npz", counts=counts, donors=np.array(donors),
                    types=np.array(types), genes=genes)
print(dt.groupby(["dataset", "diagnosis"]).size(), dt.sex.value_counts(dropna=False), flush=True)
