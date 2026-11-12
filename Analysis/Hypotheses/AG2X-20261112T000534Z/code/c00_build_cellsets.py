#!/usr/bin/env python3
"""AG2X Stage C, step C0: build the conditioning cell sets and the Stage C variant list. No hosted call.

Donor sets = cohorts (GSE256398, GSE244832, GSE202379), so the three sets are donor-disjoint by construction.
Eligible donors: the CCL1 rule (MASLD, >= 10 CCL20+ hepatocytes, >= 50 CCL20- hepatocytes), read from
CCL1 out_ccl1_job23189976/detectability_by_donor.tsv and re-derived here from the counts (assert equal).
Per set:
  ccl20pos  up to 256 CCL20+ hepatocyte nuclei (>= 1 CCL20 UMI), allocated across the set's donors by equal
            share (round robin, capped at what a donor has), random within donor.
  otherhep  the same number per donor of CCL20- hepatocyte nuclei, each the nearest unused CCL20- hepatocyte
            of the same donor on (within-donor percentile of log total UMI, of log genes detected), the
            CCL1 / DIG1 job 53 depth features; greedy in seeded random order, without replacement.
  stellate  up to 256 stellate/fibroblast nuclei from the same donors (GSE256398 label stellate_fibroblast;
            atlas label Fibroblasts), equal-share allocation, random within donor.
Counts: GSE256398 CellBender integer counts (the files CCL1 used); atlas raw X (cellranger counts). Gene IDs
are Ensembl IDs without version; duplicated IDs are summed. Seed 20261111.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
HYP = ROOT / "Analysis/Hypotheses"
PKG = HYP / "AG2X-20261112T000534Z"
C = PKG / "stageC"
CCL1 = HYP / "CCL1-20261110T235204Z/out_ccl1_job23189976"
SN_DIR = HYP / "V3-20261008T015239Z/downloads_untrusted"
SN_LAB = HYP / "DIG1-20261024T113951Z/out_snrna_labels/lineage_per_nucleus.tsv.gz"
ATLAS = ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
SEED = 20261111
N_MAX = 256
CCL20 = "ENSG00000115009"
ENS = re.compile(r"^ENS[A-Z]*G[0-9]{11}$")
SETS = {"set1": "GSE256398_snrna", "set2": "GSE244832", "set3": "GSE202379"}
rng = np.random.default_rng(SEED)


def log(*a):
    print(*a, flush=True)


def eligible():
    d = pd.read_csv(CCL1 / "detectability_by_donor.tsv", sep="\t")
    d = d[(d.disease == 1) & (d.n_ccl20_pos_hep >= 10) & (d.n_hep - d.n_ccl20_pos_hep >= 50)]
    return {c: d[d.cohort == c].set_index("donor") for c in SETS.values()}


def clean_ids(ids):
    ids = np.array([str(g).split(".")[0] for g in ids])
    bad = [g for g in ids if not ENS.match(g)]
    assert not bad, f"non-Ensembl IDs: {bad[:5]}"
    return ids


def collapse(X, ids):
    """Sum columns with the same Ensembl ID. Returns X (csr int32), unique ids, n duplicated ids."""
    u, inv = np.unique(ids, return_inverse=True)
    ndup = int(len(ids) - len(u))
    if ndup:
        M = sparse.csr_matrix((np.ones(len(ids)), (np.arange(len(ids)), inv)), shape=(len(ids), len(u)))
        X = (X @ M).tocsr()
    else:
        order = np.argsort(ids)
        X = X[:, order].tocsr()
        u = ids[order]
    X.eliminate_zeros()
    assert np.allclose(X.data, np.round(X.data)), "non-integer counts"
    return sparse.csr_matrix(X.astype(np.int32)), u, ndup


def pct(x):
    from scipy.stats import rankdata
    return (rankdata(x) - 0.5) / len(x)


def alloc(avail: dict, n=N_MAX):
    """Equal-share allocation, round robin in sorted donor order, capped at availability."""
    out = {d: 0 for d in sorted(avail)}
    left = min(n, sum(avail.values()))
    while left:
        for d in sorted(avail):
            if left and out[d] < avail[d]:
                out[d] += 1
                left -= 1
    return out


def sample_sets(meta: pd.DataFrame):
    """meta: one row per nucleus with donor, lineage (hep/stellate/other), umi, ngene, ccl20, row (int)."""
    picks = {"ccl20pos": [], "otherhep": [], "stellate": []}
    hep = meta[meta.lineage == "hep"]
    pos_av = hep[hep.ccl20 > 0].groupby("donor").size().to_dict()
    a = alloc(pos_av)
    for d in sorted(a):
        H = hep[hep.donor == d].copy()
        H["f1"] = pct(np.log(H.umi.to_numpy()))
        H["f2"] = pct(np.log(H.ngene.to_numpy()))
        P = H[H.ccl20 > 0]
        N = H[H.ccl20 == 0]
        sel = np.sort(rng.choice(len(P), size=a[d], replace=False))
        Ps = P.iloc[sel]
        tree = cKDTree(N[["f1", "f2"]].to_numpy())
        used = set()
        order = rng.permutation(len(Ps))
        match = [None] * len(Ps)
        for i in order:
            k = min(len(N), 64)
            while True:
                _, nn = tree.query(Ps[["f1", "f2"]].to_numpy()[i], k=k)
                nn = np.atleast_1d(nn)
                free = [j for j in nn if j not in used]
                if free or k >= len(N):
                    break
                k = min(len(N), k * 4)
            j = free[0]
            used.add(j)
            match[i] = j
        picks["ccl20pos"].append(Ps)
        picks["otherhep"].append(N.iloc[match])
    donors = sorted(a)
    st = meta[(meta.lineage == "stellate") & meta.donor.isin(donors)]
    sa = alloc(st.groupby("donor").size().to_dict())
    for d in sorted(sa):
        S = st[st.donor == d]
        picks["stellate"].append(S.iloc[np.sort(rng.choice(len(S), size=sa[d], replace=False))])
    return {k: pd.concat(v) if v else pd.DataFrame() for k, v in picks.items()}, a, sa


def write_set(name, state, X, ids, obs, cohort):
    import anndata as ad
    obs = obs.reset_index(drop=True)
    o = pd.DataFrame({"cohort": cohort, "donor": obs.donor.astype(str).to_numpy(), "barcode": obs.barcode.astype(str).to_numpy(),
                      "total_umi": X.sum(axis=1).A1, "n_genes": np.diff(X.indptr), "ccl20_umi": obs.ccl20.to_numpy()},
                     index=[f"c{i}" for i in range(X.shape[0])])
    assert X.shape[0] <= N_MAX and X.shape[0] > 0 and (np.diff(X.indptr) > 0).all()
    A = ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=pd.Index(ids, name="gene_id")))
    p = C / "inputs/cellsets" / f"{name}_{state}.h5ad"
    A.write_h5ad(p)
    return o


def snrna(el):
    import scanpy as sc
    labs = pd.read_csv(SN_LAB, sep="\t", index_col=0)
    parts, mats, ids0 = [], [], None
    for d in sorted(el.index):
        f = sorted(SN_DIR.glob(f"GSM*_{d}_CB_raw_feature_bc_matrix_filtered.h5"))
        assert len(f) == 1, (d, f)
        a = sc.read_10x_h5(f[0], gex_only=True)
        sym = np.asarray(a.var_names, dtype=str)
        ids = clean_ids(a.var["gene_ids"].to_numpy())
        assert sym[ids == CCL20].tolist() == ["CCL20"]
        X = sparse.csr_matrix(a.X)
        bc = np.array([f"{d}_{b}" for b in a.obs_names])
        umi = X.sum(axis=1).A1
        ng = np.diff(X.indptr)
        mt = np.char.startswith(sym, "MT-")
        pmt = X[:, mt].sum(axis=1).A1 / np.maximum(umi, 1) * 100
        keep = (ng >= 200) & (pmt < 5) & np.isin(bc, labs.index)
        X = X[keep]
        bc = bc[keep]
        X, u, ndup = collapse(X, ids)
        if ids0 is None:
            ids0 = u
        assert np.array_equal(u, ids0), "gene IDs differ between GSE256398 files"
        lin = labs.loc[bc, "lineage"].to_numpy()
        lin = np.where(lin == "hepatocyte", "hep", np.where(lin == "stellate_fibroblast", "stellate", "other"))
        ci = int(np.searchsorted(u, CCL20))
        parts.append(pd.DataFrame({"donor": d, "barcode": bc, "lineage": lin, "umi": X.sum(axis=1).A1, "ngene": np.diff(X.indptr),
                                   "ccl20": X[:, ci].toarray().ravel()}))
        mats.append(X)
        log("snrna", d, X.shape, "dup ids", ndup)
    meta = pd.concat(parts, ignore_index=True)
    meta["row"] = np.arange(len(meta))
    return meta, sparse.vstack(mats).tocsr(), ids0


def cat(o, k):
    g = o[k]
    c = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
    return c[g["codes"][:]]


def atlas(els):
    mapping = {}
    for ds in ["GSE202379", "GSE244832"]:
        t = pd.read_csv(ROOT / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
        for r in t.itertuples(index=False):
            for s in str(r.rna_srrs).split(";"):
                if s.strip():
                    mapping[s.strip()] = f"{ds}_{r.donor_id}"
    with h5py.File(ATLAS, "r") as f:
        o = f["obs"]
        ds_all = cat(o, "dataset")
        sm_all = cat(o, "sample")
        ct_all = cat(o, "cell_type")
        bc_all = np.array([x.decode() for x in o["_index"][:]])
        ids_raw = clean_ids([x.decode() for x in f["raw"]["var"]["gene_ids"][:]])
        donor_all = np.array([mapping.get(s, "") for s in sm_all])
        want = set().union(*[set(e.index) for e in els.values()])
        rows = np.where(np.isin(ds_all, ["GSE202379", "GSE244832"]) & np.isin(donor_all, list(want)))[0]
        grp = f["raw"]["X"]
        indptr = grp["indptr"][:]
        ncol = len(ids_raw)
        ci_cols = np.where(ids_raw == CCL20)[0]
        assert len(ci_cols) == 1
        umi = np.zeros(len(rows))
        ng = np.zeros(len(rows))
        cc = np.zeros(len(rows))
        pos = 0
        for b0 in range(0, len(indptr) - 1, 20_000):
            b1 = min(b0 + 20_000, len(indptr) - 1)
            sel = rows[(rows >= b0) & (rows < b1)]
            if not len(sel):
                continue
            lo, hi = indptr[b0], indptr[b1]
            B = sparse.csr_matrix((grp["data"][lo:hi], grp["indices"][lo:hi], indptr[b0:b1 + 1] - lo), shape=(b1 - b0, ncol))[sel - b0]
            k = len(sel)
            umi[pos:pos + k] = B.sum(axis=1).A1
            ng[pos:pos + k] = np.diff(B.indptr)
            cc[pos:pos + k] = B[:, ci_cols[0]].toarray().ravel()
            pos += k
        assert pos == len(rows)
        ct = ct_all[rows]
        meta = pd.DataFrame({"dataset": ds_all[rows], "donor": donor_all[rows], "barcode": bc_all[rows],
                             "lineage": np.where(ct == "Hepatocytes", "hep", np.where(ct == "Fibroblasts", "stellate", "other")),
                             "umi": umi, "ngene": ng, "ccl20": cc, "row": rows})
    return meta, ids_raw


def extract_atlas_rows(rows, ids_raw):
    """Full count rows for the selected atlas rows (in the order given)."""
    order = np.argsort(rows)
    srt = rows[order]
    out = [None] * len(rows)
    with h5py.File(ATLAS, "r") as f:
        grp = f["raw"]["X"]
        indptr = grp["indptr"]
        for i, r in zip(order, srt):
            lo, hi = indptr[r], indptr[r + 1]
            out[i] = sparse.csr_matrix((grp["data"][lo:hi], grp["indices"][lo:hi], [0, hi - lo]), shape=(1, len(ids_raw)))
    return sparse.vstack(out).tocsr()


def check_vs_ccl1(meta, el, tag):
    h = meta[meta.lineage == "hep"].groupby("donor").agg(n_hep=("ccl20", "size"), n_pos=("ccl20", lambda v: int((v > 0).sum())))
    for d, r in el.iterrows():
        assert int(h.loc[d, "n_hep"]) == int(r.n_hep) and int(h.loc[d, "n_pos"]) == int(r.n_ccl20_pos_hep), \
            f"{tag} {d}: rebuilt {h.loc[d].tolist()} vs CCL1 {[r.n_hep, r.n_ccl20_pos_hep]}"
    log(tag, "hepatocyte and CCL20+ counts equal to CCL1 for", len(el), "donors")


def variants():
    v = pd.read_csv(PKG / "stageB/inputs/variants_scored_final.tsv", sep="\t", dtype={"tier": int})
    v = v[v.tier == 1].fillna("")
    mem = v[v.roles.str.split(";").apply(lambda r: "member" in r)]
    bg = v[v.roles.str.split(";").apply(lambda r: "background" in r) & ~v.uid.isin(mem.uid)]
    r = np.random.default_rng(SEED)
    bgs = bg.iloc[np.sort(r.choice(len(bg), size=N_BG, replace=False))]
    out = pd.concat([mem.assign(stageC_role="member"), bgs.assign(stageC_role="background")]).sort_values("uid")
    pilot = sorted(r.choice(sorted(mem.uid), size=10, replace=False).tolist())
    out["pilot"] = out.uid.isin(pilot)
    out.to_csv(C / "inputs/variants_stageC.tsv", sep="\t", index=False)
    return {"tier1_total": int(len(v)), "tier1_members": int(len(mem)), "tier1_background": int(len(bg)),
            "stageC_members": int(len(mem)), "stageC_background": int(len(bgs)), "excluded_lead_proxy_control": int(len(v) - len(mem) - len(bg)),
            "pilot": pilot}


N_BG = 800


def main():
    el = eligible()
    summary = {"seed": SEED, "n_max": N_MAX, "sets": {}}
    rows_manifest = []
    # set 1: GSE256398
    meta, X, ids = snrna(el["GSE256398_snrna"])
    check_vs_ccl1(meta, el["GSE256398_snrna"], "GSE256398")
    sets = {"set1": (meta, lambda r: X[r.row.to_numpy()], ids, "GSE256398")}
    # sets 2, 3: atlas
    am, ids_raw = atlas({k: el[k] for k in ["GSE244832", "GSE202379"]})
    for s, ds in [("set2", "GSE244832"), ("set3", "GSE202379")]:
        m = am[am.dataset == ds].reset_index(drop=True)
        check_vs_ccl1(m, el[ds], ds)
        sets[s] = (m, None, ids_raw, ds)
    for s, (m, getter, gids, cohort) in sets.items():
        picks, a, sa = sample_sets(m)
        summary["sets"][s] = {"cohort": cohort, "donors": len(a)}
        for state, P in picks.items():
            if getter is not None:
                Xs = getter(P)
                Xs, u = Xs, gids
            else:
                Xs = extract_atlas_rows(P.row.to_numpy(), gids)
                Xs, u, ndup = collapse(Xs, gids)
                summary["sets"][s]["atlas_duplicated_ids"] = ndup
            o = write_set(s, state, Xs, u, P, cohort)
            for d, g in o.groupby("donor"):
                rows_manifest.append({"set": s, "cohort": cohort, "state": state, "donor": d, "n_nuclei": len(g),
                                      "median_umi": float(g.total_umi.median()), "median_genes": float(g.n_genes.median()),
                                      "ccl20_umi": float(g.ccl20_umi.sum())})
            summary["sets"][s][state] = {"n_nuclei": int(len(o)), "n_donors": int(o.donor.nunique()), "n_genes_var": int(len(u)),
                                         "median_umi": float(o.total_umi.median()), "median_genes": float(o.n_genes.median()),
                                         "nnz_per_cell_median": float(o.n_genes.median()), "ccl20_umi_total": float(o.ccl20_umi.sum())}
            log(s, state, summary["sets"][s][state])
        avail = m[m.lineage == "hep"].assign(pos=lambda x: x.ccl20 > 0).groupby("donor").pos.agg(["sum", "size"])
        for d in a:
            rows_manifest.append({"set": s, "cohort": cohort, "state": "available", "donor": d,
                                  "n_ccl20pos_hep": int(avail.loc[d, "sum"]), "n_hep": int(avail.loc[d, "size"]),
                                  "n_stellate": int(((m.donor == d) & (m.lineage == "stellate")).sum())})
    pd.DataFrame(rows_manifest).to_csv(C / "inputs/cellset_manifest_local.tsv", sep="\t", index=False)
    summary["variants"] = variants()
    json.dump(summary, open(C / "inputs/cellset_summary.json", "w"), indent=1)
    log(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
