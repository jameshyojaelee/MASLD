#!/usr/bin/env python3
"""A3 step 4: motif similarity clusters from the 879 PFMs only (no outcome data).

Similarity = max over strand and offset (overlap >= max(5, min(La, Lb) - 2) columns, capped
at min(La, Lb)) of the Pearson correlation of the flattened overlapping probability
matrices. Average linkage on 1 - similarity, cut at 0.25 (PRESPEC fallback rule if the
cluster count is outside 100-400).
"""
import os
from collections import Counter

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

A3 = os.environ.get("A3_OUT", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption")
WORK = f"{A3}/work"
RES = f"{A3}/results"
BASES = list("ACGT")


def probs():
    pf = pd.read_csv(f"{WORK}/motif_pfm.tsv", sep="\t")
    meta = pd.read_csv(f"{WORK}/motif_meta.tsv", sep="\t")
    mats = {}
    for mid, g in pf.groupby("motif_id", sort=False):
        c = g.sort_values("pos")[BASES].to_numpy(dtype=float) + 0.2
        mats[mid] = c / c.sum(axis=1, keepdims=True)
    return meta, [mats[m] for m in meta.motif_id]


def best_corr(a, b):
    La, Lb = a.shape[0], b.shape[0]
    need = min(max(5, min(La, Lb) - 2), min(La, Lb))
    best = -1.0
    for bb in (b, b[::-1, ::-1]):
        for off in range(-(Lb - need), La - need + 1):
            a0, b0 = max(0, off), max(0, -off)
            ov = min(La - a0, Lb - b0)
            if ov < need:
                continue
            x = a[a0:a0 + ov].ravel()
            y = bb[b0:b0 + ov].ravel()
            x = x - x.mean()
            y = y - y.mean()
            den = np.sqrt((x * x).sum() * (y * y).sum())
            if den > 0:
                r = (x * y).sum() / den
                if r > best:
                    best = r
    return best


def main():
    meta, P = probs()
    n = len(P)
    sim = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            sim[i, j] = sim[j, i] = best_corr(P[i], P[j])
    np.save(f"{WORK}/motif_similarity.npy", sim)
    dist = np.clip(1 - sim, 0, None)
    np.fill_diagonal(dist, 0)
    Z = linkage(squareform(dist, checks=False), method="average")
    cut = 0.25
    lab = fcluster(Z, t=cut, criterion="distance")
    k = lab.max()
    rule = "primary_cut_0.25"
    if k < 100 or k > 400:
        rule = "fallback"
        best = None
        for c in np.arange(0.10, 0.4001, 0.01):
            kk = fcluster(Z, t=c, criterion="distance").max()
            if 150 <= kk <= 300 and (best is None or abs(c - 0.25) < abs(best - 0.25)):
                best = c
        cut = float(best)
        lab = fcluster(Z, t=cut, criterion="distance")
        k = lab.max()
    print(f"cut {cut:.2f} ({rule}): {k} clusters", flush=True)

    meta["cluster_num"] = lab
    rows = []
    for c, g in meta.groupby("cluster_num"):
        fam = Counter(x.split(";")[0] for x in g.family.fillna("") if x)
        names = sorted(g.name.tolist())
        label = "/".join(names[:3]) + ("/..." if len(names) > 3 else "")
        rows.append({"cluster_num": c, "cluster_size": len(g), "cluster_label": label,
                     "cluster_family": fam.most_common(1)[0][0] if fam else "",
                     "cluster_members": ",".join(names)})
    cl = pd.DataFrame(rows)
    cl = cl.sort_values(["cluster_size", "cluster_label"], ascending=[False, True]).reset_index(drop=True)
    cl["cluster_id"] = [f"C{i + 1:03d}" for i in range(len(cl))]
    meta = meta.merge(cl[["cluster_num", "cluster_id", "cluster_label", "cluster_family", "cluster_size"]], on="cluster_num")
    meta.drop(columns="cluster_num").to_csv(f"{RES}/motif_clusters.tsv", sep="\t", index=False)
    cl.drop(columns="cluster_num").to_csv(f"{RES}/motif_cluster_summary.tsv", sep="\t", index=False)
    pd.DataFrame([{"cut_distance": cut, "rule": rule, "n_clusters": k, "n_motifs": n,
                   "n_singletons": int((cl.cluster_size == 1).sum())}]).to_csv(
        f"{RES}/motif_clustering_params.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
