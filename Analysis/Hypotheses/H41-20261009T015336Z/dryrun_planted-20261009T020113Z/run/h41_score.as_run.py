"""H41 step B: score HEPG2 shRNA knockdown consensus signatures against the refit disease signature.

Usage: python3 -I h41_score.py <H41_dir> <extract_dir> <deg_results.csv> <out_dir>

Implements PRESPEC_H41.md: hairpin-level signatures, MODZ-style gene consensus,
WTCS (primary) and Spearman (secondary) against fixed 50/50 landmark disease
sets, cell-line label-permutation null and plate-matched null (10,000 draws
each, SeedSequence(20260923).spawn(2)), BH over the 12 candidate TFs,
positive control (HNF4A top decile), negative control (12 drawn genes),
residual test beyond HNF4A, and descriptive analyses.
Environment H41_N_NULL overrides the draw count (dry runs only).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

SEED = 20260923
N_NULL = int(os.environ.get("H41_N_NULL", "10000"))
SET_SIZE = 50
MIN_SET = 20
PADJ = 0.05
CHUNK = 500
MIN_WT = 0.01
POS = "HNF4A"
NAMED = ["HNF4A", "HNF1A", "HNF1B", "FOXA1", "FOXA2", "FOXA3", "CEBPA", "CEBPB",
         "ONECUT1", "ONECUT2", "NR1H4", "PPARA", "RXRA", "ESRRA", "ESRRG", "HNF4G", "NR5A2"]
EXPECTED_PRESENT = ["HNF4A", "HNF1A", "HNF1B", "FOXA1", "FOXA2", "FOXA3", "CEBPA", "CEBPB",
                    "NR1H4", "RXRA", "ESRRA", "HNF4G"]


# ---------------------------------------------------------------- core statistics
def rank_std(x):
    """Row-wise average ranks, centred and scaled to unit norm (Spearman = dot product)."""
    x = np.atleast_2d(np.asarray(x, dtype=np.float64))
    r = rankdata(x, axis=1)
    r -= r.mean(axis=1, keepdims=True)
    nrm = np.linalg.norm(r, axis=1, keepdims=True)
    nrm[nrm == 0] = 1.0
    return r / nrm


def modz_batch(z, r):
    """MODZ-style consensus. z, r: (B, k, G) z-scores and rank_std rows. Returns (B, G), weights, mean pair rho."""
    b, k, _ = z.shape
    if k == 1:
        return z[:, 0, :], np.ones((b, 1)), np.full(b, np.nan)
    c = np.einsum("bkg,bjg->bkj", r, r)
    off = ~np.eye(k, dtype=bool)
    pair_mean = c[:, off].reshape(b, -1).mean(axis=1)
    cc = np.where(off[None], np.clip(c, 0.0, None), 0.0)
    w = cc.sum(axis=2) / (k - 1)
    w = np.maximum(w, MIN_WT)
    w = w / w.sum(axis=1, keepdims=True)
    return np.einsum("bk,bkg->bg", w, z), w, pair_mean


def es_matrix(x, mask):
    """Weighted (p = 1) KS enrichment score of gene set `mask` for every row of x (rows = signatures)."""
    x = np.atleast_2d(np.asarray(x, dtype=np.float64))
    mask = np.asarray(mask, dtype=bool)
    n_gene = x.shape[1]
    k = int(mask.sum())
    order = np.argsort(-x, axis=1, kind="stable")
    absx = np.take_along_axis(np.abs(x), order, axis=1)
    hits = mask[order]
    w = absx * hits
    nr = w.sum(axis=1)
    safe = np.where(nr > 0, nr, 1.0)
    d = np.cumsum(w, axis=1) / safe[:, None] - np.cumsum(~hits, axis=1) / (n_gene - k)
    idx = np.argmax(np.abs(d), axis=1)
    out = d[np.arange(x.shape[0]), idx]
    return np.where(nr > 0, out, 0.0)


def wtcs(es_up, es_dn):
    return np.where(np.sign(es_up) != np.sign(es_dn), (es_up - es_dn) / 2.0, 0.0)


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    ranked = pv[o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def p_upper(null, obs):
    null = np.asarray(null)
    return (1.0 + np.sum(null >= obs)) / (len(null) + 1.0)


def p_lower(null, obs):
    null = np.asarray(null)
    return (1.0 + np.sum(null <= obs)) / (len(null) + 1.0)


def percentile_among(values, i):
    """(others below + 0.5 x tied others) / (n - 1) for element i."""
    v = np.asarray(values)
    others = np.delete(v, i)
    return (np.sum(others < v[i]) + 0.5 * np.sum(others == v[i])) / len(others)


class Universe:
    """Gene universe (indices into the 978 landmark order) with fixed up/down sets and logFC ranks."""

    def __init__(self, idx, up_mask978, dn_mask978, lfc978):
        self.idx = np.asarray(idx)
        self.up = up_mask978[self.idx]
        self.dn = dn_mask978[self.idx]
        self.lfc_r = rank_std(lfc978[self.idx])[0]

    def score(self, x):
        xu = np.atleast_2d(x)[:, self.idx]
        eu = es_matrix(xu, self.up)
        ed = es_matrix(xu, self.dn)
        rho = rank_std(xu) @ self.lfc_r
        return {"wtcs": wtcs(eu, ed), "rho": rho, "down": -ed, "es_up": eu, "es_dn": ed}


def score_chunked(u, x):
    parts = [u.score(x[a:a + CHUNK]) for a in range(0, len(x), CHUNK)]
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


# ---------------------------------------------------------------- inputs
def load_disease(deg_path, ext_dir):
    lm = pd.read_csv(os.path.join(ext_dir, "landmark_order.tsv"), sep="\t", dtype=str)
    mp = pd.read_csv(os.path.join(ext_dir, "landmark_map.tsv"), sep="\t", dtype=str)
    mp = mp.set_index("pr_gene_id").reindex(lm.pr_gene_id)
    deg = pd.read_csv(deg_path)
    deg["ensembl"] = deg.gene.str.replace(r"\.[0-9]+(_PAR_Y)?$", "", regex=True)
    deg = deg.set_index("ensembl")
    ens = mp.ensembl.to_numpy()
    mapped = pd.notna(ens) & pd.Series(ens).isin(deg.index).to_numpy()
    r_idx = np.flatnonzero(mapped)
    lfc = np.full(978, np.nan)
    padj = np.full(978, np.nan)
    lfc[r_idx] = deg.loc[ens[r_idx], "logFC"].to_numpy()
    padj[r_idx] = deg.loc[ens[r_idx], "padj"].to_numpy()
    tab = pd.DataFrame({"pr_gene_id": lm.pr_gene_id, "symbol": lm.pr_gene_symbol,
                        "ensembl": ens, "logFC": lfc, "padj": padj, "in_R": mapped})
    sig = tab[tab.in_R & (tab.padj < PADJ)]
    up_c = sig[sig.logFC > 0].sort_values("logFC", ascending=False)
    dn_c = sig[sig.logFC < 0].sort_values("logFC", ascending=True)
    sizes = {"R": int(mapped.sum()), "padj_up_available": len(up_c), "padj_dn_available": len(dn_c)}
    for name, side in (("up", up_c), ("down", dn_c)):
        if len(side) < MIN_SET:
            raise SystemExit(f"STOP S4: only {len(side)} {name} genes")
    up_idx = up_c.index[:SET_SIZE].to_numpy()
    dn_idx = dn_c.index[:SET_SIZE].to_numpy()
    tab["set"] = ""
    tab.loc[up_idx, "set"] = "up"
    tab.loc[dn_idx, "set"] = "down"
    up_mask = np.zeros(978, bool)
    dn_mask = np.zeros(978, bool)
    up_mask[up_idx] = True
    dn_mask[dn_idx] = True
    sizes.update({"up": int(up_mask.sum()), "down": int(dn_mask.sum())})
    lfc0 = np.nan_to_num(lfc)
    return tab, r_idx, up_mask, dn_mask, lfc0, sizes


def hairpin_table(meta, z):
    meta = meta.reset_index(drop=True).copy()
    meta["plate"] = meta.sig_id.str.split(":").str[0]
    hp_ids = np.array(sorted(meta.pert_id.unique()))
    pos = {h: i for i, h in enumerate(hp_ids)}
    hi = meta.pert_id.map(pos).to_numpy()
    counts = np.bincount(hi, minlength=len(hp_ids)).astype(float)
    hz = np.zeros((len(hp_ids), z.shape[1]))
    np.add.at(hz, hi, z.astype(np.float64))
    hz /= counts[:, None]
    g = meta.groupby("pert_id")
    hp = pd.DataFrame({"pert_id": hp_ids,
                       "gene": g.pert_iname.first().reindex(hp_ids).to_numpy(),
                       "n_sigs": counts.astype(int),
                       "plates": g.plate.apply(lambda s: ",".join(sorted(set(s)))).reindex(hp_ids).to_numpy()})
    if (g.pert_iname.nunique() > 1).any():
        raise SystemExit("a hairpin maps to more than one gene")
    return hp, hz


# ---------------------------------------------------------------- main
def main(h41_dir, ext_dir, deg_path, out_dir):
    os.makedirs(out_dir, exist_ok=False)
    rep = {"n_null": N_NULL, "seed": SEED, "set_size": SET_SIZE}

    tab, r_idx, up_mask, dn_mask, lfc, sizes = load_disease(deg_path, ext_dir)
    rep["disease_sizes"] = sizes
    rep["flag_lt900_mapped"] = sizes["R"] < 900
    tab.to_csv(os.path.join(out_dir, "disease_signature_landmarks.tsv"), sep="\t", index=False)
    sym978 = tab.symbol.to_numpy()
    sym_pos = {s: i for i, s in enumerate(sym978)}
    r_set = set(r_idx.tolist())

    universes = {}

    def universe_for(gene):
        """Full R, or R minus the target when the target is a mapped landmark."""
        drop = sym_pos.get(gene)
        key = drop if (drop is not None and drop in r_set) else -1
        if key not in universes:
            idx = r_idx if key == -1 else r_idx[r_idx != key]
            universes[key] = Universe(idx, up_mask, dn_mask, lfc)
        return key, universes[key]

    full_key, full_u = universe_for("__none__")

    meta = pd.read_csv(os.path.join(ext_dir, "meta_trt_sh.tsv"), sep="\t", dtype=str)
    z = np.load(os.path.join(ext_dir, "z_trt_sh.npy"))
    assert len(meta) == len(z) and z.shape[1] == 978
    hp, hz = hairpin_table(meta, z)
    del z
    hr = rank_std(hz)
    n_hp = len(hp)
    gene_hp = hp.groupby("gene").indices
    genes = np.array(sorted(gene_hp))
    k_gene = np.array([len(gene_hp[g]) for g in genes])
    rep["hepg2_trt_sh"] = {"signatures": int(len(meta)), "hairpins": int(n_hp), "genes": int(len(genes)),
                           "plates": int(meta.sig_id.str.split(":").str[0].nunique())}

    present = [g for g in NAMED if g in gene_hp and len(gene_hp[g]) >= 2]
    rep["candidates_present"] = present
    rep["candidates_absent"] = [g for g in NAMED if g not in present]
    if present != EXPECTED_PRESENT:
        rep["warning_present_differs_from_prespec"] = True
    if POS not in present:
        json.dump(rep, open(os.path.join(out_dir, "run_report.json"), "w"), indent=2)
        raise SystemExit("STOP S2: HNF4A absent")
    if len(present) < 5:
        json.dump(rep, open(os.path.join(out_dir, "run_report.json"), "w"), indent=2)
        raise SystemExit("STOP S1: fewer than 5 candidates -> untestable")
    nc = pd.read_csv(os.path.join(h41_dir, "negative_control_genes.tsv"), sep="\t", dtype={"gene": str})
    nc_genes = nc.gene.tolist()
    assert all(g in gene_hp for g in nc_genes), "negative-control gene missing"
    focus = present + nc_genes

    # ---------------- observed gene-level consensus for every HEPG2 knockdown gene
    cons = {}
    obs_rows = []
    for g in genes:
        idx = gene_hp[g]
        c, w, pm = modz_batch(hz[idx][None], hr[idx][None])
        cons[g] = c[0]
        key, u = universe_for(g)
        s = u.score(c)
        obs_rows.append({"gene": g, "k": len(idx), "target_dropped_from_universe": key != -1,
                         "universe_size": len(u.idx), "wtcs": s["wtcs"][0], "rho": s["rho"][0],
                         "down": s["down"][0], "es_up": s["es_up"][0], "es_dn": s["es_dn"][0],
                         "hairpin_pair_rho": pm[0]})
    obs = pd.DataFrame(obs_rows)
    for col in ("wtcs", "rho", "down"):
        v = obs[col].to_numpy()
        obs[f"pct_{col}"] = [percentile_among(v, i) for i in range(len(v))]
        obs[f"rank_{col}"] = rankdata(-v, method="min").astype(int)

    # ---------------- cell-line null (child 0) and plate-matched null (child 1)
    ss_cell, ss_plate = np.random.SeedSequence(SEED).spawn(2)
    rng_cell = np.random.default_rng(ss_cell)
    rng_plate = np.random.default_rng(ss_plate)
    keep_k = sorted({len(gene_hp[g]) for g in focus})
    null = {}
    null_cons = {}
    for k in sorted(set(k_gene.tolist())):
        draws = np.stack([rng_cell.choice(n_hp, size=k, replace=False) for _ in range(N_NULL)])
        cs, pms = [], []
        for a in range(0, N_NULL, CHUNK):
            d = draws[a:a + CHUNK]
            c, _, pm = modz_batch(hz[d], hr[d])
            cs.append(c)
            pms.append(pm)
        c = np.concatenate(cs)
        s = score_chunked(full_u, c)
        s["pair_rho"] = np.concatenate(pms)
        null[k] = s
        if k in keep_k:
            null_cons[k] = c
    obs["p_cell_wtcs"] = [p_upper(null[k]["wtcs"], v) for k, v in zip(obs.k, obs.wtcs)]
    obs["p_cell_rho"] = [p_upper(null[k]["rho"], v) for k, v in zip(obs.k, obs.rho)]
    obs.to_csv(os.path.join(out_dir, "gene_scores_all_hepg2.tsv"), sep="\t", index=False)
    np.savez_compressed(os.path.join(out_dir, "null_cell_full_universe.npz"),
                        **{f"k{k}_{m}": null[k][m] for k in null for m in ("wtcs", "rho", "down", "pair_rho")})
    qs = [0.5, 0.9, 0.95, 0.99]
    pd.DataFrame([{"k": k, "stat": m, **{f"q{int(q * 100)}": np.quantile(null[k][m], q) for q in qs}}
                  for k in null for m in ("wtcs", "rho", "down")]).to_csv(
        os.path.join(out_dir, "null_cell_quantiles.tsv"), sep="\t", index=False)
    calib = {"all_genes": int(len(obs)),
             "frac_p_cell_wtcs_lt_0.05": float((obs.p_cell_wtcs < 0.05).mean()),
             "frac_p_cell_rho_lt_0.05": float((obs.p_cell_rho < 0.05).mean()),
             "frac_p_cell_wtcs_lt_0.05_by_k": {str(int(k)): float((p < 0.05).mean())
                                               for k, p in obs.groupby("k").p_cell_wtcs}}

    # plate pools
    hp_plate_sets = [set(p.split(",")) for p in hp.plates]
    plate_pool = {}
    for i, ps in enumerate(hp_plate_sets):
        for p in ps:
            plate_pool.setdefault(p, []).append(i)
    plate_pool = {p: np.array(v) for p, v in plate_pool.items()}

    focus_rows = []
    hairpin_rows = []
    for g in focus:
        idx = gene_hp[g]
        k = len(idx)
        key, u = universe_for(g)
        o = obs.set_index("gene").loc[g]
        # cell-line null on this gene's universe
        nc_cell = null[k] if key == -1 else score_chunked(u, null_cons[k])
        # plate-matched null
        own_plates = [sorted(hp_plate_sets[i])[0] for i in idx]
        need = pd.Series(own_plates).value_counts().sort_index()
        draws = np.empty((N_NULL, k), dtype=np.int64)
        for b in range(N_NULL):
            chosen = []
            for p, n in need.items():
                pool = plate_pool[p]
                if chosen:
                    pool = pool[~np.isin(pool, chosen)]
                chosen.extend(rng_plate.choice(pool, size=int(n), replace=False).tolist())
            draws[b] = chosen
        cs = [modz_batch(hz[draws[a:a + CHUNK]], hr[draws[a:a + CHUNK]])[0] for a in range(0, N_NULL, CHUNK)]
        sp = score_chunked(u, np.concatenate(cs))
        hs = u.score(hz[idx])
        focus_rows.append({
            "gene": g, "role": "candidate" if g in present else "negative_control", "k": k,
            "plates": ",".join(need.index), "universe_size": len(u.idx),
            "wtcs": o.wtcs, "rho": o.rho, "down": o.down, "es_up": o.es_up, "es_dn": o.es_dn,
            "pct_wtcs": o.pct_wtcs, "rank_wtcs": int(o.rank_wtcs), "pct_rho": o.pct_rho,
            "pct_down": o.pct_down,
            "p_cell_wtcs": p_upper(nc_cell["wtcs"], o.wtcs),
            "p_cell_rho": p_upper(nc_cell["rho"], o.rho),
            "p_cell_down": p_upper(nc_cell["down"], o.down),
            "p_cell_wtcs_lower": p_lower(nc_cell["wtcs"], o.wtcs),
            "p_plate_wtcs": p_upper(sp["wtcs"], o.wtcs),
            "p_plate_rho": p_upper(sp["rho"], o.rho),
            "hairpin_pair_rho": o.hairpin_pair_rho,
            "pct_pair_rho_in_null": float(np.mean(null[k]["pair_rho"] < o.hairpin_pair_rho)),
            "n_hairpins_wtcs_pos": int(np.sum(hs["wtcs"] > 0)),
            "spearman_with_HNF4A_consensus": float(rank_std(cons[g])[0] @ rank_std(cons[POS])[0]),
        })
        for j, i in enumerate(idx):
            hairpin_rows.append({"gene": g, "pert_id": hp.pert_id[i], "plates": hp.plates[i],
                                 "n_sigs": int(hp.n_sigs[i]), "wtcs": hs["wtcs"][j], "rho": hs["rho"][j],
                                 "down": hs["down"][j]})
    foc = pd.DataFrame(focus_rows)

    # hairpin percentiles among all hairpin-level signatures (full universe)
    all_hp = score_chunked(full_u, hz)
    hpf = pd.DataFrame(hairpin_rows)
    hpf["pct_wtcs_among_hairpins"] = [float(np.mean(all_hp["wtcs"] < v) + 0.5 * np.mean(all_hp["wtcs"] == v))
                                      for v in hpf.wtcs]

    # BH families
    cand = foc[foc.role == "candidate"].copy()
    negc = foc[foc.role == "negative_control"].copy()
    cand["q_T_wtcs"] = bh(cand.p_cell_wtcs)
    cand["q_T_rho"] = bh(cand.p_cell_rho)
    negc["q_N_wtcs"] = bh(negc.p_cell_wtcs)

    hn = cand.set_index("gene").loc[POS]
    pos_pass = bool(hn.pct_wtcs >= 0.90)
    neg_pass = bool((negc.q_N_wtcs < 0.05).sum() == 0)
    cand["meets_q"] = cand.q_T_wtcs < 0.05
    cand["meets_plate"] = cand.p_plate_wtcs < 0.05
    cand["called"] = cand.meets_q & cand.meets_plate & pos_pass & neg_pass & (cand.gene != POS)
    cand["plate_dependent"] = cand.meets_q & ~cand.meets_plate

    # beyond HNF4A residual test
    beyond = []
    x = np.column_stack([np.ones(978), cons[POS]])
    for g in cand.loc[cand.called, "gene"]:
        k = len(gene_hp[g])
        key, u = universe_for(g)
        beta, *_ = np.linalg.lstsq(x, cons[g], rcond=None)
        res = cons[g] - x @ beta
        so = u.score(res)
        bn, *_ = np.linalg.lstsq(x, null_cons[k].T, rcond=None)
        rn = null_cons[k] - (x @ bn).T
        sn = score_chunked(u, rn)
        beyond.append({"gene": g, "wtcs_residual": so["wtcs"][0], "p_residual": p_upper(sn["wtcs"], so["wtcs"][0]),
                       "slope_on_HNF4A": beta[1]})
    beyond = pd.DataFrame(beyond, columns=["gene", "wtcs_residual", "p_residual", "slope_on_HNF4A"])
    if len(beyond):
        beyond["q_residual"] = bh(beyond.p_residual)
    beyond.to_csv(os.path.join(out_dir, "beyond_hnf4a.tsv"), sep="\t", index=False)

    # Broad CGS signatures (descriptive)
    cm = pd.read_csv(os.path.join(ext_dir, "meta_trt_sh.cgs.tsv"), sep="\t", dtype=str)
    cz = np.load(os.path.join(ext_dir, "z_trt_sh.cgs.npy")).astype(np.float64)
    cgs_all = score_chunked(full_u, cz)
    cgs_rows = []
    for g in present:
        ii = np.flatnonzero(cm.pert_iname.to_numpy() == g)
        for i in ii:
            key, u = universe_for(g)
            v = u.score(cz[i])["wtcs"][0]
            others = np.delete(cgs_all["wtcs"], i)
            cgs_rows.append({"gene": g, "cgs_sig_id": cm.sig_id[i], "cgs_wtcs": v,
                             "cgs_pct_wtcs": float((np.sum(others < v) + 0.5 * np.sum(others == v)) / len(others))})
    cgs = pd.DataFrame(cgs_rows, columns=["gene", "cgs_sig_id", "cgs_wtcs", "cgs_pct_wtcs"])
    cgs.to_csv(os.path.join(out_dir, "cgs_candidates.tsv"), sep="\t", index=False)

    # CEBPA on-target check (only landmark candidate)
    ontarget = []
    for g in present:
        if g in sym_pos:
            j = sym_pos[g]
            allz = hz[:, j]
            for i in gene_hp[g]:
                ontarget.append({"gene": g, "pert_id": hp.pert_id[i], "target_z": hz[i, j],
                                 "pct_target_z_among_hairpins": float(np.mean(allz < hz[i, j]))})
            ontarget.append({"gene": g, "pert_id": "consensus", "target_z": cons[g][j],
                             "pct_target_z_among_hairpins": np.nan})
    pd.DataFrame(ontarget).to_csv(os.path.join(out_dir, "on_target_landmark.tsv"), sep="\t", index=False)

    cand.to_csv(os.path.join(out_dir, "candidate_results.tsv"), sep="\t", index=False)
    negc.to_csv(os.path.join(out_dir, "negative_control_results.tsv"), sep="\t", index=False)
    hpf.to_csv(os.path.join(out_dir, "hairpin_scores_focus.tsv"), sep="\t", index=False)

    called = cand.loc[cand.called, "gene"].tolist()
    if not (pos_pass and neg_pass):
        answer = "not interpretable in this assay"
    elif called:
        answer = "yes"
    else:
        answer = "no"
    summary = {
        "answer": answer,
        "called": called,
        "plate_dependent": cand.loc[cand.plate_dependent, "gene"].tolist(),
        "family_T_size": int(len(cand)),
        "family_N_size": int(len(negc)),
        "positive_control": {"gene": POS, "pct_wtcs": float(hn.pct_wtcs), "rank_wtcs": int(hn.rank_wtcs),
                             "n_genes": int(len(obs)), "p_cell_wtcs": float(hn.p_cell_wtcs),
                             "threshold": 0.90, "pass": pos_pass},
        "negative_control": {"n_q_lt_0.05": int((negc.q_N_wtcs < 0.05).sum()),
                             "n_p_lt_0.05": int((negc.p_cell_wtcs < 0.05).sum()), "pass": neg_pass},
        "calibration": calib,
    }
    rep.update(summary)
    rep["universes"] = {str(k): len(v.idx) for k, v in universes.items()}
    with open(os.path.join(out_dir, "run_report.json"), "w") as fh:
        json.dump(rep, fh, indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
