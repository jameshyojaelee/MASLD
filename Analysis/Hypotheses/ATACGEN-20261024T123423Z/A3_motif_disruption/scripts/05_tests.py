#!/usr/bin/env python3
"""A3 step 5: Q1 enrichment (stratified permutation), Q2 risk-allele direction (loci as
unit), Q3 convergence with chromVAR disease shifts, Q4 top examples. See PRESPEC.md."""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

A3 = os.environ.get("A3_OUT", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption")
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WORK = f"{A3}/work"
RES = f"{A3}/results"
CHROMVAR = f"{BASE}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/chromvar/chromvar_results.tsv.gz"
R2 = f"{BASE}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
LINEAGES = ["hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"]
B = 10_000
B_LOLO = 2_000
GAP = 250_000
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    o = np.argsort(pv)
    r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty_like(pv)
    qq[o] = np.minimum(r, 1)
    q[ok] = qq
    return q


# ---------------------------------------------------------------- permutation engine
def perm_test(strata, w, D, n_perm, seed, chunk=250):
    """Stratified permutation of weights w (0/1 labels or PIPs) within strata.
    Returns observed, analytic E and Var, and two-sided Monte Carlo p for each column of D."""
    order = np.argsort(strata, kind="stable")
    s = strata[order]
    ws = w[order].astype(np.float32)
    Ds = D[order].astype(np.float32)
    obs = ws @ Ds
    # analytic moments of the permutation distribution
    uniq, start, cnt = np.unique(s, return_index=True, return_counts=True)
    E = np.zeros(D.shape[1])
    V = np.zeros(D.shape[1])
    for st, n in zip(start, cnt):
        wb = ws[st:st + n].astype(float)
        db = Ds[st:st + n].astype(float)
        E += wb.mean() * db.sum(0)
        if n > 1:
            V += ((wb - wb.mean()) ** 2).sum() * ((db - db.mean(0)) ** 2).sum(0) / (n - 1)
    rng = np.random.default_rng(seed)
    ge = np.zeros(D.shape[1])
    le = np.zeros(D.shape[1])
    tol = 1e-4 * np.maximum(1, np.abs(obs))
    sf = s.astype(np.float64)
    done = 0
    while done < n_perm:
        k = min(chunk, n_perm - done)
        key = sf[None, :] + rng.random((k, s.size))
        perm = np.argsort(key, axis=1)
        Wp = np.empty((k, s.size), dtype=np.float32)
        np.put_along_axis(Wp, perm, np.broadcast_to(ws, (k, s.size)), axis=1)
        T = Wp @ Ds
        ge += (T >= obs - tol).sum(0)
        le += (T <= obs + tol).sum(0)
        done += k
    p = np.minimum(1, 2 * np.minimum((1 + ge) / (n_perm + 1), (1 + le) / (n_perm + 1)))
    return obs.astype(float), E, V, p


def informative(strata, is_case, is_ctrl):
    df = pd.DataFrame({"s": strata, "a": is_case, "b": is_ctrl})
    g = df.groupby("s").agg(a=("a", "sum"), b=("b", "sum"))
    good = g.index[(g.a > 0) & (g.b > 0)]
    return np.isin(strata, good)


# ---------------------------------------------------------------- data
def load():
    # keep_default_na=False: pandas would otherwise read the vclass label "null" as NaN;
    # R fwrite writes missing values as empty fields
    v = pd.read_csv(f"{WORK}/variants_hg38.tsv.gz", sep="\t", dtype={"key": np.int64}, low_memory=False,
                    keep_default_na=False, na_values=[""])
    assert v.vclass.isin(["high", "member_low", "null"]).all()
    keys = pd.read_csv(f"{WORK}/scored_keys.tsv.gz", sep="\t", dtype={"key": np.int64}).key.to_numpy()
    vp = v[v.in_any_peak].reset_index(drop=True)
    assert np.array_equal(vp.key.to_numpy(), keys)
    s_ref = np.load(f"{WORK}/s_ref.npy")
    s_alt = np.load(f"{WORK}/s_alt.npy")
    meta = pd.read_csv(f"{WORK}/motif_meta_thresholds.tsv", sep="\t")
    cl = pd.read_csv(f"{RES}/motif_clusters.tsv", sep="\t")
    meta = meta.merge(cl[["motif_id", "cluster_id", "cluster_label", "cluster_family", "cluster_size"]], on="motif_id")
    return v, vp, s_ref, s_alt, meta


def locus_clusters(v):
    """Single linkage of gated-member hg38 positions (gap <= GAP); Null variants join the
    cluster of their nearest gated member within GAP."""
    mem = v[v.vclass.isin(["high", "member_low"])].sort_values(["chr", "pos38"])
    out = pd.Series(index=v.index, dtype=object)
    for c, g in mem.groupby("chr"):
        pos = g.pos38.to_numpy()
        cid = np.cumsum(np.r_[True, np.diff(pos) > GAP])
        first = {k: pos[cid == k].min() for k in np.unique(cid)}
        last = {k: pos[cid == k].max() for k in np.unique(cid)}
        nm = np.array([f"chr{c}:{first[k]}-{last[k]}" for k in cid], dtype=object)
        out[g.index] = nm
        sel = (v.chr == c) & (v.vclass == "null")
        q = v.pos38[sel].to_numpy()
        j = np.searchsorted(pos, q)
        jl = np.clip(j - 1, 0, len(pos) - 1)
        jr = np.clip(j, 0, len(pos) - 1)
        dl = np.abs(q - pos[jl])
        dr = np.abs(pos[jr] - q)
        best = np.where(dl <= dr, jl, jr)
        out[v.index[sel]] = np.where(np.minimum(dl, dr) <= GAP, nm[best], None)
    return out


def tss_bin(d):
    d = np.asarray(d, dtype=float)
    return np.where(d < 2000, 0, np.where(d < 20000, 1, 2))


def qbin(x, k, ref):
    """Quantile bins 0..k-1 from cut points of ref; missing values get bin k."""
    x = np.asarray(x, dtype=float)
    cuts = np.nanquantile(np.asarray(ref, dtype=float), np.linspace(0, 1, k + 1)[1:-1])
    return np.where(np.isnan(x), k, np.digitize(x, cuts))


def maf_bin(m):
    m = np.asarray(m, dtype=float)
    return np.where(np.isnan(m), 3, np.where(m < 0.05, 0, np.where(m < 0.2, 1, 2)))


def factor(*cols):
    if len(cols[0]) == 0:
        return np.zeros(0, dtype=np.int64)
    df = pd.DataFrame({f"c{i}": c for i, c in enumerate(cols)})
    return pd.factorize(pd.MultiIndex.from_frame(df.astype(str)))[0]


# ---------------------------------------------------------------- main
def main():
    v, vp, s_ref, s_alt, meta = load()
    motif_ids = meta.motif_id.tolist()
    thr = meta["threshold_p5e-5"].to_numpy()
    clusters = meta.drop_duplicates("cluster_id").sort_values("cluster_id")[["cluster_id", "cluster_label", "cluster_family", "cluster_size"]].reset_index(drop=True)
    cid = clusters.cluster_id.tolist()
    M = np.zeros((len(motif_ids), len(cid)), dtype=bool)
    cidx = {c: j for j, c in enumerate(cid)}
    for i, c in enumerate(meta.cluster_id):
        M[i, cidx[c]] = True

    v["locus"] = locus_clusters(v)
    vp = vp.merge(v[["key", "locus"]], on="key", how="left")
    hit = np.maximum(s_ref, s_alt) >= thr[None, :]
    delta = s_alt - s_ref
    absd = np.abs(delta)
    dis = {d: hit & (absd >= d) for d in (1.0, 2.0, 3.0)}
    cdis = {d: (dis[d].astype(np.float32) @ M.astype(np.float32)) > 0 for d in dis}

    # cluster delta for strong disruptions: member with the largest |delta| among strong members
    strong = dis[3.0]
    cl_delta = np.full((len(vp), len(cid)), np.nan, dtype=np.float32)
    cl_best_motif = np.full((len(vp), len(cid)), -1, dtype=np.int32)
    for j in range(len(cid)):
        mi = np.nonzero(M[:, j])[0]
        a = np.where(strong[:, mi], absd[:, mi], -1)
        b = a.argmax(axis=1)
        okr = a.max(axis=1) >= 0
        cl_delta[okr, j] = delta[okr, mi[b[okr]]]
        cl_best_motif[okr, j] = mi[b[okr]]

    summary = {"n_variants_in_any_peak": int(len(vp)), "n_motifs": len(motif_ids), "n_clusters": len(cid)}
    counts = []
    q1_rows, q1w_rows, glob_rows, sens_rows, lolo_rows = [], [], [], [], []
    q1_by_lineage = {}

    for li, L in enumerate(LINEAGES, start=1):
        inL = vp[f"peak_{L}"].to_numpy(bool)
        has_locus = vp.locus.notna().to_numpy()
        high = (vp.vclass == "high").to_numpy()
        null = (vp.vclass == "null").to_numpy()
        memb = vp.vclass.isin(["high", "member_low"]).to_numpy()
        base = inL & has_locus & (high | null)
        idx = np.nonzero(base)[0]
        gc_t = qbin(vp.gc101.to_numpy()[idx], 3, vp.gc101.to_numpy()[idx])
        tb = tss_bin(vp.tss_dist.to_numpy()[idx])
        st = factor(vp.locus.to_numpy()[idx], gc_t, tb)
        inf = informative(st, high[idx], null[idx])
        ii = idx[inf]
        st_i = factor(st[inf])
        nh, nn = int(high[ii].sum()), int(null[ii].sum())
        counts.append({"lineage": L, "high_in_peak": int((inL & high).sum()), "null_in_peak_with_locus": int((inL & null & has_locus).sum()),
                       "high_informative": nh, "null_informative": nn, "n_strata_informative": int(st_i.max() + 1) if len(st_i) else 0,
                       "n_locus_clusters_informative": len(set(vp.locus.to_numpy()[ii][high[ii]]))})
        if nh < 20:
            q1_by_lineage[L] = None
            glob_rows.append({"lineage": L, "analysis": "binary", "testable": False})
            continue
        Dc = cdis[2.0][ii].astype(np.float32)
        Dm = dis[2.0][ii].astype(np.float32)
        Dt = Dc.sum(1, keepdims=True)
        D = np.hstack([Dc, Dm, Dt])
        w = high[ii].astype(np.float32)
        O, E, Var, p = perm_test(st_i, w, D, B, 20261024 + li)
        pooled = D.sum(0)
        nc, nm = len(cid), len(motif_ids)
        loc_i = vp.locus.to_numpy()[ii]

        def contrib(col):
            sel = (w > 0) & (D[:, col] > 0)
            return sorted(set(loc_i[sel]))

        tab = pd.DataFrame({
            "lineage": L, "level": ["cluster"] * nc + ["motif"] * nm,
            "id": cid + motif_ids,
            "label": clusters.cluster_label.tolist() + meta.name.tolist(),
            "family": clusters.cluster_family.tolist() + meta.family.fillna("").tolist(),
            "cluster_id": cid + meta.cluster_id.tolist(),
            "n_high_informative": nh, "n_null_informative": nn,
            "pooled_disrupted": pooled[:-1], "O": O[:-1], "E": E[:-1], "Var": Var[:-1], "p_mc": p[:-1],
        })
        tab["O_over_E"] = tab.O / tab.E.replace(0, np.nan)
        tab["z"] = (tab.O - tab.E) / np.sqrt(tab.Var.replace(0, np.nan))
        tab["testable"] = tab.pooled_disrupted >= 10
        tab["q_bh"] = np.nan
        for lev in ("cluster", "motif"):
            s = (tab.level == lev) & tab.testable
            tab.loc[s, "q_bh"] = bh(tab.loc[s, "p_mc"].to_numpy())
        tab["n_loci_contributing"] = [len(contrib(j)) for j in range(nc + nm)]
        tab["direction"] = np.where(tab.O_over_E > 1, "enriched", np.where(tab.O_over_E < 1, "depleted", "equal"))
        q1_rows.append(tab)
        q1_by_lineage[L] = tab[tab.level == "cluster"].set_index("id")
        glob_rows.append({"lineage": L, "analysis": "binary", "testable": True, "n_high": nh, "n_null": nn,
                          "O": O[-1], "E": E[-1], "O_over_E": O[-1] / E[-1], "z": (O[-1] - E[-1]) / np.sqrt(Var[-1]),
                          "p_mc": p[-1], "mean_clusters_disrupted_high": O[-1] / nh,
                          "mean_clusters_disrupted_null": float(Dt[w == 0].mean())})

        # PIP-weighted (gated members vs Null)
        basew = inL & has_locus & (memb | null)
        idw = np.nonzero(basew)[0]
        gcw = qbin(vp.gc101.to_numpy()[idw], 3, vp.gc101.to_numpy()[idw])
        stw = factor(vp.locus.to_numpy()[idw], gcw, tss_bin(vp.tss_dist.to_numpy()[idw]))
        infw = informative(stw, memb[idw], null[idw])
        iw = idw[infw]
        ww = np.where(memb[iw], vp.max_gated_pip.to_numpy()[iw], 0).astype(np.float32)
        Dw = np.hstack([cdis[2.0][iw].astype(np.float32), cdis[2.0][iw].sum(1, keepdims=True).astype(np.float32)])
        Ow, Ew, Vw, pw = perm_test(factor(stw[infw]), ww, Dw, B, 20261024 + 10 + li)
        tw = pd.DataFrame({"lineage": L, "id": cid, "label": clusters.cluster_label, "family": clusters.cluster_family,
                           "n_members_informative": int(memb[iw].sum()), "sum_pip": float(ww.sum()),
                           "pooled_disrupted": Dw[:, :-1].sum(0), "O_w": Ow[:-1], "E_w": Ew[:-1], "Var_w": Vw[:-1], "p_mc": pw[:-1]})
        tw["O_over_E"] = tw.O_w / tw.E_w.replace(0, np.nan)
        tw["z"] = (tw.O_w - tw.E_w) / np.sqrt(tw.Var_w.replace(0, np.nan))
        tw["testable"] = tw.pooled_disrupted >= 10
        tw.loc[tw.testable, "q_bh"] = bh(tw.loc[tw.testable, "p_mc"].to_numpy())
        q1w_rows.append(tw)
        glob_rows.append({"lineage": L, "analysis": "pip_weighted", "testable": True, "n_high": int(memb[iw].sum()),
                          "n_null": int(null[iw].sum()), "O": Ow[-1], "E": Ew[-1], "O_over_E": Ow[-1] / Ew[-1],
                          "z": (Ow[-1] - Ew[-1]) / np.sqrt(Vw[-1]), "p_mc": pw[-1]})

        # sensitivities for clusters with primary q < 0.10, plus the global statistic
        sel_c = tab[(tab.level == "cluster") & (tab.q_bh < 0.10)].id.tolist()
        cols = [cidx[c] for c in sel_c]
        if True:
            def run(stv, wv, Dmat, seed):
                o, e, va, pp = perm_test(stv, wv, Dmat, B, seed)
                return o, e, pp
            sens_specs = []
            # B: genome-wide strata (no locus): all in-peak High + Null in fine-mapped windows
            ib = np.nonzero(inL & (high | null))[0]
            stb = factor(qbin(vp.gc101.to_numpy()[ib], 5, vp.gc101.to_numpy()[ib]), tss_bin(vp.tss_dist.to_numpy()[ib]),
                         maf_bin(vp.maf.to_numpy()[ib]), qbin(vp[f"signal_{L}"].to_numpy()[ib], 3, vp[f"signal_{L}"].to_numpy()[ib]))
            infb = informative(stb, high[ib], null[ib])
            sens_specs.append(("B_genomewide_gc5_tss_maf_signal", ib[infb], factor(stb[infb]), high[ib[infb]], cdis[2.0]))
            # C: EUR-only High
            high_eur = (vp.vclass == "high").to_numpy() & (vp.max_gated_pip_eur.fillna(0).to_numpy() >= 0.1)
            ic0 = np.nonzero(inL & has_locus & (high_eur | null))[0]
            stc = factor(vp.locus.to_numpy()[ic0], qbin(vp.gc101.to_numpy()[ic0], 3, vp.gc101.to_numpy()[ic0]), tss_bin(vp.tss_dist.to_numpy()[ic0]))
            infc = informative(stc, high_eur[ic0], null[ic0])
            if high_eur[ic0[infc]].sum() >= 20:
                sens_specs.append(("C_eur_only_high", ic0[infc], factor(stc[infc]), high_eur[ic0[infc]], cdis[2.0]))
            # D: delta thresholds
            sens_specs.append(("D_delta1", ii, st_i, high[ii], cdis[1.0]))
            sens_specs.append(("D_delta3", ii, st_i, high[ii], cdis[3.0]))
            # E: locus-only strata
            ste = factor(vp.locus.to_numpy()[idx])
            infe = informative(ste, high[idx], null[idx])
            sens_specs.append(("E_locus_only", idx[infe], factor(ste[infe]), high[idx[infe]], cdis[2.0]))
            for k, (name, rows, stv, wv, Dsrc) in enumerate(sens_specs):
                Dmat = Dsrc[rows][:, cols].astype(np.float32) if cols else np.zeros((len(rows), 0), np.float32)
                Dmat = np.hstack([Dmat, Dsrc[rows].sum(1, keepdims=True).astype(np.float32)])
                o, e, pp = run(stv, wv.astype(np.float32), Dmat, 20261024 + 100 * (k + 1) + li)
                for j, c in enumerate(sel_c + ["GLOBAL"]):
                    sens_rows.append({"lineage": L, "id": c, "sensitivity": name, "n_high": int(wv.sum()), "O": o[j], "E": e[j],
                                      "O_over_E": o[j] / e[j] if e[j] > 0 else np.nan, "p_mc": pp[j]})
            # F: leave-one-locus-out for selected clusters
            for c in sel_c:
                j = cidx[c]
                for lo in contrib(j):
                    keep = loc_i != lo
                    stk = factor(st_i[keep])
                    infk = informative(stk, w[keep] > 0, w[keep] == 0)
                    Dk = Dc[keep][infk][:, [j]]
                    o, e, va, pp = perm_test(factor(stk[infk]), w[keep][infk], Dk, B_LOLO, 20261024 + 1_000_000 + 10_000 * li + j)
                    lolo_rows.append({"lineage": L, "id": c, "dropped_locus": lo, "O": o[0], "E": e[0],
                                      "O_over_E": o[0] / e[0] if e[0] > 0 else np.nan, "p_mc": pp[0]})

    q1 = pd.concat(q1_rows, ignore_index=True)
    q1w = pd.concat(q1w_rows, ignore_index=True)
    glob = pd.DataFrame(glob_rows)
    for a in ("binary", "pip_weighted"):
        s = (glob.analysis == a) & glob.testable.astype(bool)
        glob.loc[s, "q_bh_across_lineages"] = bh(glob.loc[s, "p_mc"].to_numpy())
    sens = pd.DataFrame(sens_rows)
    lolo = pd.DataFrame(lolo_rows)
    pd.DataFrame(counts).to_csv(f"{RES}/q1_informative_counts.tsv", sep="\t", index=False)
    q1[q1.level == "cluster"].drop(columns="cluster_id").to_csv(f"{RES}/q1_cluster_enrichment.tsv", sep="\t", index=False)
    q1[q1.level == "motif"].to_csv(f"{RES}/q1_motif_enrichment.tsv", sep="\t", index=False)
    q1w.to_csv(f"{RES}/q1_weighted_cluster_enrichment.tsv", sep="\t", index=False)
    glob.to_csv(f"{RES}/q1_global.tsv", sep="\t", index=False)
    sens.to_csv(f"{RES}/q1_sensitivity.tsv", sep="\t", index=False)
    lolo.to_csv(f"{RES}/q1_leave_one_locus_out.tsv", sep="\t", index=False)

    # robust flag
    qc = q1[(q1.level == "cluster")].copy()
    rob = []
    for _, r in qc[qc.q_bh < 0.05].iterrows():
        s = sens[(sens.lineage == r.lineage) & (sens.id == r.id)].set_index("sensitivity")
        same = lambda nm: (nm in s.index) and (s.loc[nm, "p_mc"] < 0.05) and (np.sign(s.loc[nm, "O_over_E"] - 1) == np.sign(r.O_over_E - 1))
        lo = lolo[(lolo.lineage == r.lineage) & (lolo.id == r.id)]
        rob.append({"lineage": r.lineage, "id": r.id, "label": r.label, "O_over_E": r.O_over_E, "q_bh": r.q_bh,
                    "B_same_dir_p05": same("B_genomewide_gc5_tss_maf_signal"), "E_same_dir_p05": same("E_locus_only"),
                    "lolo_max_p": lo.p_mc.max() if len(lo) else np.nan, "lolo_min_OE": lo.O_over_E.min() if len(lo) else np.nan})
    rob = pd.DataFrame(rob)
    if len(rob):
        rob["robust"] = rob.B_same_dir_p05 & rob.E_same_dir_p05 & (rob.lolo_max_p < 0.05)
    rob.to_csv(f"{RES}/q1_robustness.tsv", sep="\t", index=False)

    # ------------------------------------------------------------ Q2
    risk = vp.risk_allele.astype(object).to_numpy()
    flip = vp.strand_flipped.to_numpy(bool)
    risk = np.array([COMP.get(a, a) if f and isinstance(a, str) else a for a, f in zip(risk, flip)], dtype=object)
    rsign = np.where(risk == vp.alt.to_numpy(), 1.0, np.where(risk == vp.ref.to_numpy(), -1.0, np.nan))
    pip = vp.max_gated_pip.to_numpy()
    high = (vp.vclass == "high").to_numpy()
    nullv = (vp.vclass == "null").to_numpy() & vp.locus.notna().to_numpy()
    pal = vp.palindromic.to_numpy(bool)
    loc = vp.locus.to_numpy()

    def pairs(sel, use_risk=True):
        r, c = np.nonzero(np.isfinite(cl_delta) & sel[:, None])
        d = cl_delta[r, c].astype(float)
        if use_risk:
            d = d * rsign[r]
        ok = np.isfinite(d)
        return pd.DataFrame({"row": r[ok], "cluster": np.array(cid)[c[ok]], "d": d[ok], "locus": loc[r[ok]],
                             "w": np.where(np.isfinite(pip[r[ok]]), pip[r[ok]], 1.0)})

    def vote(df, by=("locus",)):
        if len(df) == 0:
            return dict(n_loci=0, n_weaken=0, n_strengthen=0, frac_weaken=np.nan, p_binom=np.nan)
        g = (df.w * np.sign(df.d)).groupby([df[b] for b in by]).sum()
        g = g[g != 0]
        k = int((g < 0).sum())
        n = int(len(g))
        return dict(n_loci=n, n_weaken=k, n_strengthen=n - k, frac_weaken=k / n if n else np.nan,
                    p_binom=stats.binomtest(k, n, 0.5).pvalue if n else np.nan,
                    n_pairs=int(len(df)), n_variants=int(df.row.nunique()),
                    variant_frac_weaken=float((df.groupby("row").apply(lambda x: (x.w * np.sign(x.d)).sum()) < 0).mean()))

    q2 = []
    anyp = vp.in_any_peak.to_numpy(bool)
    q2.append({"set": "pooled_all_lineages", "orientation": "risk", **vote(pairs(high & anyp))})
    for L in LINEAGES:
        q2.append({"set": L, "orientation": "risk", **vote(pairs(high & vp[f"peak_{L}"].to_numpy(bool)))})
    q2 = pd.DataFrame(q2)
    s = q2.set != "pooled_all_lineages"
    q2.loc[s, "q_bh_across_lineages"] = bh(q2.loc[s, "p_binom"].to_numpy())
    ctrl = [
        {"control": "risk_is_alt", **vote(pairs(high & anyp & (rsign == 1)))},
        {"control": "risk_is_ref", **vote(pairs(high & anyp & (rsign == -1)))},
        {"control": "high_alt_orientation", **vote(pairs(high & anyp, use_risk=False))},
        {"control": "null_alt_orientation", **vote(pairs(nullv & anyp, use_risk=False).assign(w=1.0))},
        {"control": "risk_excluding_palindromic", **vote(pairs(high & anyp & ~pal))},
    ]
    ctrl = pd.DataFrame(ctrl)
    hv = vp[high & anyp]
    multi = hv[hv.n_gated_rows >= 2]
    ctrl_conc = {"n_high_in_peak": int(len(hv)), "frac_risk_is_alt": float(np.nanmean(rsign[high & anyp] == 1)),
                 "n_multi_row": int(len(multi)), "frac_single_risk_allele": float((multi.n_distinct_risk_alleles == 1).mean()) if len(multi) else np.nan}
    q2.to_csv(f"{RES}/q2_direction_global.tsv", sep="\t", index=False)
    ctrl.to_csv(f"{RES}/q2_controls.tsv", sep="\t", index=False)
    pd.DataFrame([ctrl_conc]).to_csv(f"{RES}/q2_risk_allele_concordance.tsv", sep="\t", index=False)
    pc = pairs(high & anyp)
    rows = []
    for c, g in pc.groupby("cluster"):
        r = vote(g)
        rows.append({"cluster_id": c, **r})
    q2c = pd.DataFrame(rows).merge(clusters, on="cluster_id", how="left")
    q2c["testable"] = q2c.n_loci >= 8
    q2c.loc[q2c.testable, "q_bh"] = bh(q2c.loc[q2c.testable, "p_binom"].to_numpy())
    q2c.sort_values("p_binom").to_csv(f"{RES}/q2_direction_cluster.tsv", sep="\t", index=False)
    pairs_out = pc.assign(variant=vp.key.to_numpy()[pc.row], chr=vp.chr.to_numpy()[pc.row], pos38=vp.pos38.to_numpy()[pc.row])
    pairs_out.drop(columns="row").to_csv(f"{RES}/q2_strong_disruption_pairs.tsv.gz", sep="\t", index=False)

    # ------------------------------------------------------------ Q3
    cv = pd.read_csv(CHROMVAR, sep="\t")
    cv = cv[cv.contrast_testable.astype(str).str.upper() == "TRUE"].dropna(subset=["pvalue"])
    cv["z"] = np.sign(cv.effect) * stats.norm.isf(cv.pvalue / 2)
    cvm = meta[["motif_id", "cluster_id"]]
    q3_rows, q3_cl = [], []
    size_strata = clusters.set_index("cluster_id").cluster_size.map(lambda n: 0 if n == 1 else 1 if n <= 3 else 2 if n <= 7 else 3)
    for li, L in enumerate(LINEAGES, start=1):
        c = cv[cv.lineage == L]
        if q1_by_lineage.get(L) is None or len(c) == 0:
            q3_rows.append({"lineage": L, "testable": False, "reason": "no chromVAR contrast" if len(c) == 0 else "Q1 untestable"})
            continue
        cohorts = sorted(c.cohort.unique())
        zm = c.groupby("motif_id").z.agg(lambda x: x.sum() / np.sqrt(len(x))).rename("Z")
        qmin = c.groupby("motif_id").qvalue.min().rename("qmin")
        mm = cvm.merge(zm, on="motif_id").merge(qmin, on="motif_id")
        mm["absZ"] = mm.Z.abs()
        best = mm.sort_values("absZ", ascending=False).drop_duplicates("cluster_id").set_index("cluster_id")
        dset = mm.groupby("cluster_id").qmin.min() < 0.05
        t1 = q1_by_lineage[L]
        t1 = t1[t1.testable]
        ids = [i for i in t1.index if i in best.index]
        zq1 = t1.loc[ids, "z"].to_numpy(float)
        sabs = best.loc[ids, "absZ"].to_numpy(float)
        ok = np.isfinite(zq1) & np.isfinite(sabs)
        ids = list(np.array(ids)[ok])
        zq1, sabs = zq1[ok], sabs[ok]
        rz = stats.rankdata(zq1)
        rs = stats.rankdata(sabs)
        rho = np.corrcoef(rz, rs)[0, 1]
        strata = size_strata.loc[ids].to_numpy()
        rng = np.random.default_rng(20261124 + li)
        cnt = 0
        for _ in range(B):
            rp = rs.copy()
            for sv in np.unique(strata):
                m = strata == sv
                rp[m] = rng.permutation(rp[m])
            cnt += np.corrcoef(rz, rp)[0, 1] >= rho - 1e-12
        p_rho = (1 + cnt) / (B + 1)
        Dl = dset.reindex(ids).fillna(False).to_numpy(bool)
        qv = t1.loc[ids, "q_bh"].to_numpy(float)
        oe = t1.loc[ids, "O_over_E"].to_numpy(float)
        pv = t1.loc[ids, "p_mc"].to_numpy(float)
        R = (qv < 0.10) & (oe > 1)
        r_def = "q<0.10"
        if R.sum() < 3:
            R = (pv < 0.05) & (oe > 1)
            r_def = "exploratory_nominal_p<0.05"
        ov, p_ov = np.nan, np.nan
        if Dl.sum() >= 3 and R.sum() >= 1:
            ov = int((R & Dl).sum())
            cnt = 0
            for _ in range(B):
                dp = Dl.copy()
                for sv in np.unique(strata):
                    m = strata == sv
                    dp[m] = rng.permutation(dp[m])
                cnt += (R & dp).sum() >= ov
            p_ov = (1 + cnt) / (B + 1)
        q3_rows.append({"lineage": L, "testable": True, "chromvar_cohorts": ",".join(cohorts), "n_clusters": len(ids),
                        "spearman_rho": rho, "p_perm_one_sided": p_rho, "n_D": int(Dl.sum()), "R_definition": r_def,
                        "n_R": int(R.sum()), "overlap_R_D": ov, "p_overlap_one_sided": p_ov,
                        "overlap_expected": float(R.sum() * Dl.mean()) if len(Dl) else np.nan})
        for i, cidv in enumerate(ids):
            b = best.loc[cidv]
            q3_cl.append({"lineage": L, "cluster_id": cidv, "q1_z": zq1[i], "q1_O_over_E": oe[i], "q1_q": qv[i],
                          "chromvar_best_motif": b.motif_id, "chromvar_Z": b.Z, "chromvar_absZ": b.absZ,
                          "chromvar_min_q_member": float(mm[mm.cluster_id == cidv].qmin.min()), "in_D": bool(Dl[i]), "in_R": bool(R[i])})
    q3 = pd.DataFrame(q3_rows)
    s = q3.testable.astype(bool)
    q3.loc[s, "q_bh_across_lineages"] = bh(q3.loc[s, "p_perm_one_sided"].to_numpy())
    q3.to_csv(f"{RES}/q3_convergence.tsv", sep="\t", index=False)
    q3c = pd.DataFrame(q3_cl).merge(clusters, on="cluster_id", how="left")
    # descriptive: Q2 risk direction for clusters in both R and D
    q3c = q3c.merge(q2c[["cluster_id", "n_loci", "frac_weaken"]].rename(columns={"n_loci": "q2_n_loci", "frac_weaken": "q2_frac_weaken"}),
                    on="cluster_id", how="left")
    q3c.to_csv(f"{RES}/q3_cluster_shift.tsv", sep="\t", index=False)

    # ------------------------------------------------------------ Q4
    hi = np.nonzero(high & anyp)[0]
    sc = np.where(np.isfinite(cl_delta[hi]), np.abs(cl_delta[hi]), 0) * pip[hi][:, None]
    bestc = sc.argmax(axis=1)
    has = sc.max(axis=1) > 0
    ex = []
    for k in np.nonzero(has)[0]:
        r = hi[k]
        j = bestc[k]
        mi = cl_best_motif[r, j]
        lins = [L for L in LINEAGES if vp.loc[r, f"peak_{L}"]]
        d = float(cl_delta[r, j])
        ex.append({"chr": int(vp.chr[r]), "pos38": int(vp.pos38[r]), "pos19": int(vp.pos19[r]), "ref": vp.ref[r], "alt": vp.alt[r],
                   "max_gated_pip": pip[r], "assigned_study": vp.assigned_study[r], "trait": vp.assigned_trait[r],
                   "risk_allele_hg38": risk[r], "locus_cluster": loc[r], "lineages_with_peak": ",".join(lins),
                   "cluster_id": cid[j], "cluster_label": clusters.cluster_label[j], "best_motif": motif_ids[mi],
                   "best_motif_name": meta.name[mi], "s_ref": float(s_ref[r, mi]), "s_alt": float(s_alt[r, mi]),
                   "delta_alt_minus_ref": d, "threshold": float(thr[mi]),
                   "risk_effect": ("weakens" if d * rsign[r] < 0 else "strengthens") if np.isfinite(rsign[r]) else "risk_allele_unresolved",
                   "score_pip_x_absdelta": float(sc[k, j]), "nearest_tss_gene": vp.nearest_tss_gene[r], "tss_dist": vp.tss_dist[r],
                   "n_clusters_strong": int(np.isfinite(cl_delta[r]).sum())})
    ex = pd.DataFrame(ex).sort_values("score_pip_x_absdelta", ascending=False).reset_index(drop=True)
    sp = pd.read_csv(f"{R2}/susie_coloc_signal_pairs.tsv.gz", sep="\t")
    sp = sp[(sp.eligible.astype(str).str.upper() == "TRUE") & (sp["PP.H4"] >= 0.5) & sp.hit1.str.match(r"^\d+:\d+$")]
    sp[["hchr", "hpos"]] = sp.hit1.str.split(":", expand=True).astype(int)
    ab = pd.read_csv(f"{R2}/susie_coloc_all_gwas.csv", usecols=["gwas_name", "gene", "PP.H4.abf", "top_snp"])
    ab = ab[(ab["PP.H4.abf"] >= 0.5) & ab.top_snp.astype(str).str.match(r"^\d+:\d+")]
    ab[["tchr", "tpos"]] = ab.top_snp.str.split(":", expand=True).iloc[:, :2].astype(int)

    def genes_near(c, p):
        g1 = sp[(sp.hchr == c) & (np.abs(sp.hpos - p) <= GAP)]
        g1 = g1.sort_values("PP.H4", ascending=False).drop_duplicates("gene")
        g2 = ab[(ab.tchr == c) & (np.abs(ab.tpos - p) <= GAP)]
        g2 = g2.sort_values("PP.H4.abf", ascending=False).drop_duplicates("gene")
        s1 = ";".join(f"{a}({b:.2f},{s})" for a, b, s in zip(g1.gene, g1["PP.H4"], g1.gwas))
        s2 = ";".join(f"{a}({b:.2f})" for a, b in zip(g2.gene, g2["PP.H4.abf"]))
        lead = ";".join(sorted(set(g1[(g1.hpos == p)].gene)))
        return s1, s2, lead

    top = ex.head(30).copy()
    gg = [genes_near(c, p) for c, p in zip(top.chr, top.pos19)]
    top["coloc_r2_susie_genes_250kb"] = [g[0] for g in gg]
    top["coloc_r2_abf_genes_250kb"] = [g[1] for g in gg]
    top["coloc_r2_susie_genes_with_variant_as_gwas_lead"] = [g[2] for g in gg]
    # chromVAR and Q1 context for the example's cluster in its lineages
    ctx = []
    for _, r in top.iterrows():
        parts = []
        for L in r.lineages_with_peak.split(","):
            t = q1_by_lineage.get(L)
            if t is not None and r.cluster_id in t.index:
                parts.append(f"{L}:O/E={t.loc[r.cluster_id, 'O_over_E']:.2f},q={t.loc[r.cluster_id, 'q_bh']:.2g}")
        ctx.append(";".join(parts))
    top["q1_cluster_context"] = ctx
    cvb = cv[cv.motif_id.isin(top.best_motif)]
    top["chromvar_best_motif_per_cohort"] = [";".join(f"{a}/{b}:eff={e:.2f},q={q:.2g}" for a, b, e, q in
                                                      zip(g.cohort, g.lineage, g.effect, g.qvalue))
                                             for g in [cvb[(cvb.motif_id == m) & (cvb.lineage.isin(l.split(",")))] for m, l in zip(top.best_motif, top.lineages_with_peak)]]
    top.to_csv(f"{RES}/q4_top_examples.tsv", sep="\t", index=False)
    ex.to_csv(f"{RES}/q4_all_strong_examples.tsv.gz", sep="\t", index=False)

    summary.update({"counts": counts})
    with open(f"{RES}/summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    print("done", flush=True)


if __name__ == "__main__":
    main()
