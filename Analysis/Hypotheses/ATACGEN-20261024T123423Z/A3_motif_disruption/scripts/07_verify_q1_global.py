#!/usr/bin/env python3
"""A3 check: recompute the Q1 global O, E (exact stratified expectation) and O/E per lineage
from the per-variant table and allele scores, written independently of 05_tests.py.
Also recomputes BH for the primary global family and the hepatocyte cluster family.
Writes results/verify_q1_global_recompute.tsv. No permutation (p values are not recomputed)."""
import numpy as np
import pandas as pd

A3 = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption"
W, R = f"{A3}/work", f"{A3}/results"
LINEAGES = ["hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"]

v = pd.read_csv(f"{W}/variants_hg38.tsv.gz", sep="\t", keep_default_na=False, na_values=[""], low_memory=False,
                usecols=["key", "vclass", "chr", "pos38", "gc101", "tss_dist", "in_any_peak"] + [f"peak_{L}" for L in LINEAGES])
v["chr"] = v["chr"].astype(str)

# locus clusters: single linkage of gated members (high, member_low), gap <= 250 kb
loc = pd.Series(None, index=v.index, dtype=object)
for c, g in v[v.vclass.isin(["high", "member_low"])].groupby("chr"):
    g = g.sort_values("pos38")
    p = g.pos38.to_numpy()
    lab = np.cumsum(np.r_[1, np.diff(p) > 250_000])
    loc[g.index] = [f"{c}_{k}" for k in lab]
    nl = v[(v.chr == c) & (v.vclass == "null") & (v.in_any_peak.astype(str).str.upper() == "TRUE")]
    for i, q in zip(nl.index, nl.pos38.to_numpy()):
        d = np.abs(p - q)
        j = d.argmin()  # ties go to the lower position, as in 05
        if d[j] <= 250_000:
            loc[i] = f"{c}_{lab[j]}"
v["locus"] = loc

vp = v[v.in_any_peak.astype(str).str.upper() == "TRUE"].reset_index(drop=True)
keys = pd.read_csv(f"{W}/scored_keys.tsv.gz", sep="\t").key.to_numpy()
assert np.array_equal(vp.key.to_numpy(), keys)
s_ref = np.load(f"{W}/s_ref.npy", mmap_mode="r")
s_alt = np.load(f"{W}/s_alt.npy", mmap_mode="r")
meta = pd.read_csv(f"{W}/motif_meta_thresholds.tsv", sep="\t")
cl = pd.read_csv(f"{R}/motif_clusters.tsv", sep="\t").set_index("motif_id").loc[meta.motif_id]
thr = meta["threshold_p5e-5"].to_numpy()
cids = sorted(cl.cluster_id.unique())

# number of clusters disrupted per variant (|delta| >= 2 and either allele >= threshold)
n_clus = np.zeros(len(vp), dtype=np.int32)
clus_dis = np.zeros((len(vp), len(cids)), dtype=bool)
for a in range(0, len(vp), 5000):
    sr, sa = np.asarray(s_ref[a:a + 5000]), np.asarray(s_alt[a:a + 5000])
    d = (np.maximum(sr, sa) >= thr) & (np.abs(sa - sr) >= 2)
    for j, c in enumerate(cids):
        clus_dis[a:a + 5000, j] = d[:, (cl.cluster_id == c).to_numpy()].any(1)
n_clus = clus_dis.sum(1)

rows = []
for L in LINEAGES:
    s = vp[(vp[f"peak_{L}"].astype(str).str.upper() == "TRUE") & vp.locus.notna() & vp.vclass.isin(["high", "null"])].copy()
    s["nclus"] = n_clus[s.index]
    cuts = np.quantile(s.gc101, [1 / 3, 2 / 3])
    s["gct"] = np.digitize(s.gc101, cuts)
    s["tssb"] = np.select([s.tss_dist < 2000, s.tss_dist < 20000], [0, 1], 2)
    s["h"] = (s.vclass == "high").astype(int)
    g = s.groupby(["locus", "gct", "tssb"])
    s["nh_s"], s["n_s"] = g.h.transform("sum"), g.h.transform("size")
    inf = s[(s.nh_s > 0) & (s.nh_s < s.n_s)]
    O = float((inf.h * inf.nclus).sum())
    gi = inf.groupby(["locus", "gct", "tssb"])
    E = float((gi.h.mean() * gi.nclus.sum()).sum())
    rows.append({"lineage": L, "n_high": int(inf.h.sum()), "n_null": int((1 - inf.h).sum()), "n_strata": gi.ngroups,
                 "O": O, "E": E, "O_over_E": O / E})
out = pd.DataFrame(rows)

glob = pd.read_csv(f"{R}/q1_global.tsv", sep="\t").query("analysis == 'binary'").set_index("lineage")
out = out.join(glob[["n_high", "n_null", "O", "E", "O_over_E", "q_bh_across_lineages"]].add_suffix("_05"), on="lineage")


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty_like(p)
    q[o] = np.minimum(r, 1)
    return q


out["q_bh_recomputed_from_p05"] = bh(glob.loc[out.lineage, "p_mc"].to_numpy())
ce = pd.read_csv(f"{R}/q1_cluster_enrichment.tsv", sep="\t").query("lineage == 'hepatocyte' and testable")
out["hep_cluster_family_n"] = len(ce)
out["hep_cluster_max_abs_q_diff"] = float(np.abs(bh(ce.p_mc.to_numpy()) - ce.q_bh.to_numpy()).max())
out.to_csv(f"{R}/verify_q1_global_recompute.tsv", sep="\t", index=False)
print(out.to_string())
