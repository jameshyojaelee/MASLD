#!/usr/bin/env python3
"""Stage A statistics exactly as fixed in PRESPEC_AG2X.md (sha256 a6943f7d...; sections 2-5).

No model, calibration or prior is fitted on AlphaGenome outputs: only rank, sign and pairwise statistics.
Unit of resampling: 1-Mb block (chr~floor(pos_hg38/1e6)); BlockBootstrap of c1e2_05_recompute_v2,
10,000 draws, seed 20261111; every statistic of one analysis is computed in the same draws.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

sys.path.insert(0, str(L.REPO / "scripts/analysis/alphagenome_program"))
sys.path.insert(0, str(L.REPO / "scripts/analysis/alphagenome_atlas"))
import c1e2_05_recompute_v2 as C  # noqa: E402
import lib_atlas as la  # noqa: E402

N_BOOT = 10_000
SEED = L.SEED
RAW, TAB, IN = L.PKG / "raw", L.PKG / "tables", L.PKG / "inputs"
C1 = L.REPO / "GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z/tables"
TRUTH = L.REPO / "GWAS/finemapping/results/seqfunc/broadaway_benchmark_truth.tsv"
READOUTS = {"BL_RNA": ("bulk_liver", "RNA_SEQ"), "HEP_SCRNA": ("hepatocyte", "SC_GENE_EXPRESSION"),
            "BL_ATAC": ("bulk_liver", "ATAC"), "BL_DNASE": ("bulk_liver", "DNASE"), "HEP_SCATAC": ("hepatocyte", "SC_ATAC")}
EXPR, ACC = ["BL_RNA", "HEP_SCRNA"], ["BL_ATAC", "BL_DNASE", "HEP_SCATAC"]
METRICS, CONTRASTS, NOTES = [], [], {}


def block(uid):
    c, p = uid.split(":")[:2]
    return f"{c}~{int(p) // 1_000_000}"


# ------------------------------------------------------------------ scores
def load_v2():
    cols = ["variant_id", "gene_id", "output_type", "track_name", "imputation_mode", "raw_score"]
    var, gene, tracks = {}, {}, []
    for cname in ("bulk_liver", "hepatocyte"):
        d = pd.concat([pd.read_parquet(p, columns=cols) for p in sorted((RAW / f"v2_{cname}").glob("chunk_*.parquet"))],
                      ignore_index=True)
        t = d.groupby(["output_type", "imputation_mode", "track_name"]).size().rename("rows").reset_index()
        t.insert(0, "condition", cname)
        tracks.append(t)
        d = d[d.imputation_mode == "observed"].copy()
        d["uid"] = d.variant_id.str.replace(">", ":", regex=False)
        v = d[d.output_type.isin(["ATAC", "DNASE", "SC_ATAC"])]
        var[cname] = v.groupby(["uid", "output_type"]).raw_score.mean().unstack()
        g = d[d.output_type.isin(["RNA_SEQ", "SC_GENE_EXPRESSION"])].copy()
        g["ensg"] = g.gene_id.str.split(".").str[0]
        gene[cname] = g.groupby(["uid", "ensg", "output_type"]).raw_score.mean().unstack()
        L.log(f"v2 {cname}: {d.uid.nunique()} variants with observed rows")
    tracks = pd.concat(tracks, ignore_index=True)
    tracks["used"] = tracks.imputation_mode.eq("observed") & [
        any(c == cn and o == ot for cn, ot in READOUTS.values()) for c, o in zip(tracks.condition, tracks.output_type)]
    tracks.to_csv(TAB / "observed_tracks_used.tsv", sep="\t", index=False)
    vs, gs = {}, {}
    for r, (cn, ot) in READOUTS.items():
        if ot in ("RNA_SEQ", "SC_GENE_EXPRESSION"):
            gs[r] = gene[cn][ot] if ot in gene[cn] else pd.Series(dtype=float)
        else:
            vs[r] = var[cn][ot] if ot in var[cn] else pd.Series(dtype=float)
    return vs, gs, tracks


def load_v1():
    d = pd.concat([pd.read_parquet(p, columns=["variant_id", "gene_id", "ontology_curie", "track_name", "raw_score"])
                   for p in sorted((RAW / "v1_rna_liver").glob("chunk_*.parquet"))], ignore_index=True)
    d["uid"] = d.variant_id.str.replace(">", ":", regex=False)
    d["ensg"] = d.gene_id.str.split(".").str[0]
    NOTES["v1_hosted_liver_rna_tracks"] = sorted(d.track_name.unique().tolist())
    return d.groupby(["uid", "ensg"]).raw_score.mean()


# ------------------------------------------------------------------ statistics
def auroc(s, y):
    ok = np.isfinite(s)
    s, y = s[ok], y[ok]
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def boot(df, fns):
    """fns: name -> f(df_sorted, idx). Returns {name: (point, draws)}."""
    df = df.assign(_b=df.block_1mb.astype(str)).sort_values("_b", kind="mergesort").reset_index(drop=True)
    bb = C.BlockBootstrap(df._b.values, N_BOOT, SEED)
    allidx = np.arange(len(df))
    out = {k: (f(df, allidx), np.full(N_BOOT, np.nan)) for k, f in fns.items()}
    for i, idx in enumerate(bb.draws()):
        for k, f in fns.items():
            out[k][1][i] = f(df, idx)
    return out, int(df._b.nunique())


def rec(endpoint, readout, metric, res, n, blocks, null, extra=None):
    pt, dr = res
    lo, hi = C.ci(dr)
    p, note = C.two_sided_boot_p(dr - null) if np.isfinite(null) else (np.nan, "")
    METRICS.append({"endpoint": endpoint, "readout": readout, "metric": metric, "n": n, "blocks": blocks,
                    "estimate": pt, "ci_lo": lo, "ci_hi": hi, "null": null, "boot_p": p, "p_note": note, **(extra or {})})


def contrast(cid, endpoint, a, b, metric, ra, rb, n, blocks, family):
    d = ra[1] - rb[1]
    lo, hi = C.ci(d)
    p, note = C.two_sided_boot_p(d)
    CONTRASTS.append({"id": cid, "family": family, "endpoint": endpoint, "contrast": f"{a} - {b}", "metric": metric,
                      "n": n, "blocks": blocks, "estimate": ra[0] - rb[0], "ci_lo": lo, "ci_hi": hi, "boot_p": p, "p_note": note})


def sc(df, col, lab):
    return lambda d, i: C.signed_spearman(d[col].values[i], d[lab].values[i])


def wsc(df, col, lab, w):
    return lambda d, i: C.weighted_spearman(d[col].values[i], d[lab].values[i], d[w].values[i])


def sagree(col, lab):
    return lambda d, i: C.sign_concordance(d[col].values[i], d[lab].values[i])


# ------------------------------------------------------------------ endpoints
def spearman_endpoint(name, df, label, weight, models, contrasts_spec):
    gone = [x for x in models if df[x].isna().all()]
    if gone:
        NOTES[f"{name}_untestable_readouts"] = gone
    models = [x for x in models if x not in gone]
    n0 = len(df)
    df = df.dropna(subset=models).copy()
    NOTES[f"{name}_rows_dropped_for_missing_score"] = n0 - len(df)
    fns = {}
    for m in models:
        fns[f"{m}|unweighted"] = sc(df, m, label)
        if weight:
            fns[f"{m}|weighted"] = wsc(df, m, label, weight)
    res, nb = boot(df, fns)
    for k, r in res.items():
        m, kind = k.split("|")
        rec(name, m, f"signed_spearman_{kind}", r, len(df), nb, 0.0)
    for cid, a, b, fam in contrasts_spec:
        if a in gone or b in gone:
            continue
        kind = "weighted" if weight else "unweighted"
        contrast(cid, name, a, b, f"signed_spearman_{kind}", res[f"{a}|{kind}"], res[f"{b}|{kind}"], len(df), nb, fam)
    return res


def a1(vs):
    d = pd.read_csv(C1 / "endpoint1_matched_tierA4.tsv.gz", sep="\t")
    d["uid"] = [f"{c}:{p}:{r}:{a}" for c, p, r, a in zip(d.chr_x, d.pos_hg38, d.ref, d.alt)]
    for r in ACC:
        d[r] = d.uid.map(vs[r])
    NOTES["A1_missing"] = {r: int(d[r].isna().sum()) for r in ACC}
    models = ACC + ["alphagenome_atac_liver", "control_allele_identity", "control_position"]
    return spearman_endpoint("A1", d, "beta_alt", None, models,
                             [("S1", "BL_ATAC", "alphagenome_atac_liver", "secondary"), ("S7", "HEP_SCATAC", "BL_ATAC", "secondary")])


def a2(vs):
    d = pd.read_csv(C1 / "endpoint2_matched.tsv", sep="\t")
    for r in ACC:
        d[r] = d.uid.map(vs[r])
    NOTES["A2_missing"] = {r: int(d[r].isna().sum()) for r in ACC}
    models = ACC + ["alphagenome_atac_liver", "control_allele_identity", "control_position"]
    return spearman_endpoint("A2_495", d, "mean_log2_alt_over_ref", "precision_weight", models,
                             [("S2", "BL_ATAC", "alphagenome_atac_liver", "secondary"), ("S8", "HEP_SCATAC", "BL_ATAC", "secondary")])


def a2hep(vs):
    d = pd.read_csv(TAB / "a2_hepatocyte_lineage_v1.tsv", sep="\t")
    for r in ACC:
        d[r] = d.uid.map(vs[r])
    NOTES["A2hep_missing"] = {r: int(d[r].isna().sum()) for r in ACC}
    full = spearman_endpoint("A2_hep", d, "mean_log2_alt_over_ref", "precision_weight",
                             ACC + ["control_allele_identity", "control_position"],
                             [("S9", "HEP_SCATAC", "BL_ATAC", "secondary")])
    spearman_endpoint("A2_hep_with_v1", d, "mean_log2_alt_over_ref", "precision_weight", ["BL_ATAC", "alphagenome_atac_liver"],
                      [("S3", "BL_ATAC", "alphagenome_atac_liver", "secondary")])
    return full


def a3a(gs):
    d = pd.read_csv(TAB / "a3a_gate_set_v1.tsv", sep="\t")
    d["block_1mb"] = d.variant_uid.map(block)
    for r in EXPR:
        d[r] = [gs[r].get((u, e), np.nan) for u, e in zip(d.variant_uid, d.ensembl)]
        d[f"o_{r}"] = d.orientation_sign * d[r]
    d["o_v1"] = d.orientation_sign * d.v1_atlas_liver_rna_ref_to_alt
    d["meas"] = np.where(d.label == 1, 1.0, -1.0)
    NOTES["A3a_missing"] = {r: int(d[r].isna().sum()) for r in EXPR}
    for r in EXPR:
        s = d.dropna(subset=[r])
        if not len(s):
            continue
        res, nb = boot(s, {"auroc": lambda x, i, c=f"o_{r}": auroc(x[c].values[i], x.label.values[i]),
                           "sign": sagree(f"o_{r}", "meas")})
        rec("A3a", r, "auroc", res["auroc"], len(s), nb, 0.5)
        rec("A3a", r, "sign_agreement", res["sign"], len(s), nb, 0.5,
            {"marginal_expected": la.marginal_expected_concordance(s[f"o_{r}"], s.meas)})
    s = d.dropna(subset=["BL_RNA", "o_v1"])
    res, nb = boot(s, {"BL_RNA": lambda x, i: auroc(x.o_BL_RNA.values[i], x.label.values[i]),
                       "v1": lambda x, i: auroc(x.o_v1.values[i], x.label.values[i])})
    rec("A3a_with_v1", "v1_atlas_liver_rna", "auroc", res["v1"], len(s), nb, 0.5)
    contrast("S4", "A3a_with_v1", "BL_RNA", "v1_atlas_liver_rna", "auroc", res["BL_RNA"], res["v1"], len(s), nb, "secondary")


def orientation_qc():
    tr = pd.read_csv(TRUTH, sep="\t", dtype=str)
    tr = tr[tr.strand_ambiguous != "TRUE"].copy()
    tr["ensg"] = tr.ensembl.str.split(".").str[0]
    tr["pos_hg19"] = tr.pos_hg19.astype(int)
    tr = tr.drop_duplicates(["ensg", "pos_hg19"]).rename(columns={"pos_hg38": "truth_pos_hg38"})
    cand = pd.read_csv(IN / "eqtl_pip_candidates_hg38.tsv.gz", sep="\t",
                       usecols=["ensg", "pos_hg19", "chr_hg38", "pos_hg38", "NEA", "EA", "Beta", "liftover_hits", "strand_hg38"])
    m = tr.merge(cand, on=["ensg", "pos_hg19"])
    m = m[(m.liftover_hits == 1) & (m.strand_hg38 == "+") & (m.NEA.str.len() == 1) & (m.EA.str.len() == 1)]
    import pyfaidx
    fa = pyfaidx.Fasta(str(L.FASTA), as_raw=True, sequence_always_upper=True)
    m["fref"] = [fa[c][int(p) - 1] for c, p in zip(m.chr_hg38, m.pos_hg38)]
    ours = np.sign(m.Beta.astype(float)) * np.where(m.fref == m.NEA, 1, np.where(m.fref == m.EA, -1, 0))
    theirs = m.eqtl_sign.astype(float) * np.where(m.effect_allele.str.upper() == m.fref, -1, 1)
    ok = (ours != 0) & (theirs != 0)
    res = {"pairs_compared": int(ok.sum()), "agreement": float((ours[ok] == theirs[ok]).mean()) if ok.sum() else None,
           "hg38_position_agrees_with_truth": float((m.pos_hg38.astype(int) == m.truth_pos_hg38.astype(int)).mean())
           if len(m) else None}
    NOTES["orientation_qc"] = res
    return res


def a3b(gs, v1):
    t = pd.read_csv(IN / "a3b_a4_tops.tsv", sep="\t")
    t = t[(t.pip >= 0.9) & ~t.palindromic].copy()
    t["meas"] = np.sign(t.Beta)
    for r in EXPR:
        t[r] = [gs[r].get((u, e), np.nan) for u, e in zip(t.uid, t.ensg)]
    t["v1"] = [v1.get((u, e), np.nan) for u, e in zip(t.uid, t.ensg)]
    for c in EXPR + ["v1"]:
        t[f"ea_{c}"] = t.ea_is_alt * t[c]
    NOTES["A3b_missing"] = {c: int(t[c].isna().sum()) for c in EXPR + ["v1"]}
    NOTES["A3b_zero_predictions"] = {c: int((t[c] == 0).sum()) for c in EXPR + ["v1"]}
    out = {}
    for c in EXPR + ["v1"]:
        s = t.dropna(subset=[c])
        if not len(s):
            continue
        res, nb = boot(s, {"sign": sagree(f"ea_{c}", "meas")})
        rec("A3b", c, "sign_agreement", res["sign"], len(s), nb, 0.5,
            {"marginal_expected": la.marginal_expected_concordance(s[f"ea_{c}"][s[f"ea_{c}"] != 0], s.meas[s[f"ea_{c}"] != 0])})
        out[c] = res["sign"]
    for cid, a, b in (("S5", "BL_RNA", "v1"), ("S10", "HEP_SCRNA", "BL_RNA")):
        s = t.dropna(subset=[a, b])
        if not len(s):
            continue
        res, nb = boot(s, {a: sagree(f"ea_{a}", "meas"), b: sagree(f"ea_{b}", "meas")})
        contrast(cid, "A3b", a, b if b != "v1" else "v1_hosted_liver_rna", "sign_agreement", res[a], res[b], len(s), nb, "secondary")
    return out


def ap(score, y):
    o = np.argsort(-score, kind="mergesort")
    y = y[o]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    return float(prec[y == 1].mean()) if y.sum() else np.nan


def a4(gs, v1):
    p = pd.read_csv(IN / "a4_pairs.tsv", sep="\t")
    sc_ = {**{r: gs[r] for r in EXPR}, "v1": v1}
    for c, s in sc_.items():
        p[f"t_{c}"] = [s.get((u, e), np.nan) for u, e in zip(p.top_uid, p.ensg)]
        p[f"c_{c}"] = [s.get((u, e), np.nan) for u, e in zip(p.ctrl_uid, p.ensg)]
        at, ac = p[f"t_{c}"].abs(), p[f"c_{c}"].abs()
        p[f"i_{c}"] = np.where(at > ac, 1.0, np.where(at == ac, 0.5, 0.0))
        p.loc[at.isna() | ac.isna(), f"i_{c}"] = np.nan
    p["i_dist"] = np.where(p.top_dist_tss < p.ctrl_dist_tss, 1.0, np.where(p.top_dist_tss == p.ctrl_dist_tss, 0.5, 0.0))
    NOTES["A4_pairs_missing"] = {c: int(p[f"i_{c}"].isna().sum()) for c in sc_}
    mean = lambda col: (lambda d, i: float(np.nanmean(d[col].values[i])))  # noqa: E731
    for c in EXPR + ["v1"]:
        s = p.dropna(subset=[f"i_{c}"]).copy()
        if not len(s):
            continue
        s["delta"] = s[f"i_{c}"] - s.i_dist
        res, nb = boot(s, {"auroc": mean(f"i_{c}"), "dist": mean("i_dist"), "delta": mean("delta")})
        ex = {"genes": int(s.ensg.nunique())}
        rec("A4", c, "within_gene_auroc", res["auroc"], len(s), nb, 0.5, ex)
        rec("A4", c, "distance_baseline_auroc", res["dist"], len(s), nb, 0.5, ex)
        rec("A4", c, "auroc_minus_distance", res["delta"], len(s), nb, 0.0, ex)
        u = pd.concat([s[["ensg", "top_uid", f"t_{c}"]].drop_duplicates(["ensg", "top_uid"]).set_axis(["ensg", "uid", "s"], axis=1).assign(y=1),
                       s[["ensg", "ctrl_uid", f"c_{c}"]].set_axis(["ensg", "uid", "s"], axis=1).assign(y=0)])
        METRICS.append({"endpoint": "A4", "readout": c, "metric": "auprc_descriptive", "n": len(u), "blocks": nb,
                        "estimate": ap(u.s.abs().values, u.y.values), "null": float(u.y.mean()), "genes": ex["genes"]})
    for cid, a, b in (("S6", "BL_RNA", "v1"), ("S11", "HEP_SCRNA", "BL_RNA")):
        s = p.dropna(subset=[f"i_{a}", f"i_{b}"])
        if not len(s):
            continue
        res, nb = boot(s, {a: mean(f"i_{a}"), b: mean(f"i_{b}")})
        contrast(cid, "A4", a, b if b != "v1" else "v1_hosted_liver_rna", "within_gene_auroc", res[a], res[b], len(s), nb, "secondary")


# ------------------------------------------------------------------ multiplicity and decisions
def holm(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    adj = np.empty_like(ps)
    run = 0.0
    for k, j in enumerate(o):
        run = max(run, min(1.0, (len(ps) - k) * ps[j]))
        adj[j] = run
    return adj


def bh(ps):
    ps = np.asarray(ps, float)
    n = len(ps)
    o = np.argsort(ps)
    adj = np.empty(n)
    prev = 1.0
    for k in range(n - 1, -1, -1):
        j = o[k]
        prev = min(prev, ps[j] * n / (k + 1))
        adj[j] = prev
    return adj


def main():
    vs, gs, tracks = load_v2()
    v1 = load_v1()
    res_a1 = a1(vs)
    a2(vs)
    a2hep(vs)
    a3a(gs)
    qc = orientation_qc()
    L.log(f"orientation QC {qc}")
    if qc["agreement"] is not None and qc["agreement"] >= 0.99:
        a3b(gs, v1)
        a4(gs, v1)
    else:
        NOTES["A3b_A4_stopped"] = "orientation QC below 99% or no shared pairs"
    m = pd.DataFrame(METRICS)
    c = pd.DataFrame(CONTRASTS)

    def get(ep, r, met):
        s = m[(m.endpoint == ep) & (m.readout == r) & (m.metric == met)]
        return s.iloc[0] if len(s) else None

    dec, holm_rows = [], []
    for r in EXPR:
        p1, p2 = get("A3b", r, "sign_agreement"), get("A4", r, "auroc_minus_distance")
        if p1 is None or p2 is None:
            dec.append({"readout": r, "kind": "expression", "pass": None, "reason": "not testable or stopped"})
            continue
        adj = holm([p1.boot_p, p2.boot_p])
        holm_rows += [("A3b", r, "sign_agreement", adj[0]), ("A4", r, "auroc_minus_distance", adj[1])]
        ok = bool(p1.ci_lo > 0.5 and p2.ci_lo > 0.0)
        dec.append({"readout": r, "kind": "expression", "test1": "A3b sign agreement CI above 0.5", "test1_ci_lo": p1.ci_lo,
                    "test2": "A4 AUROC minus distance CI above 0", "test2_ci_lo": p2.ci_lo, "pass": ok,
                    "p1": p1.boot_p, "p2": p2.boot_p})
    for r in ACC:
        t1, t2 = get("A1", r, "signed_spearman_unweighted"), get("A2_hep", r, "signed_spearman_weighted")
        if t1 is None or t2 is None:
            dec.append({"readout": r, "kind": "accessibility", "pass": None, "reason": "not testable (no observed track)"})
            continue
        ok = bool(t1.ci_lo > 0 and t2.ci_lo > 0)
        dec.append({"readout": r, "kind": "accessibility", "test1": "A1 Spearman CI above 0", "test1_ci_lo": t1.ci_lo,
                    "test2": "A2-hep weighted Spearman CI above 0", "test2_ci_lo": t2.ci_lo, "pass": ok,
                    "p1": t1.boot_p, "p2": t2.boot_p})
    d = pd.DataFrame(dec)
    tested = d[d.p1.notna()] if "p1" in d else d.iloc[0:0]
    allp = np.r_[tested.p1.values, tested.p2.values].astype(float)
    if len(allp):
        adj = holm(allp)
        k = len(tested)
        d.loc[tested.index, "holm_all_tests_both_survive"] = (adj[:k] < 0.05) & (adj[k:] < 0.05)
    m["holm_within_readout"] = np.nan
    for ep, r, met, a in holm_rows:
        m.loc[(m.endpoint == ep) & (m.readout == r) & (m.metric == met), "holm_within_readout"] = a
    if len(c):
        sec = c.family == "secondary"
        c.loc[sec, "bh_q"] = bh(c.loc[sec, "boot_p"].values)
    m.to_csv(TAB / "stageA_metrics.tsv", sep="\t", index=False)
    c.to_csv(TAB / "stageA_contrasts.tsv", sep="\t", index=False)
    d.to_csv(TAB / "stageA_decisions.tsv", sep="\t", index=False)
    json.dump({"notes": NOTES, "n_boot": N_BOOT, "seed": SEED, "utc": L.utc()}, open(TAB / "stageA_notes.json", "w"),
              indent=1, default=str)
    pd.set_option("display.width", 250)
    print(m.to_string())
    print(c.to_string())
    print(d.to_string())
    print(json.dumps(NOTES, indent=1, default=str))


if __name__ == "__main__":
    main()
