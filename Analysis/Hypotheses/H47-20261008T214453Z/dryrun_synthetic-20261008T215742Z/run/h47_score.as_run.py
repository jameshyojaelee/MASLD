"""H47 step B: disease-signature reversal in LINCS L1000 liver lines, class calibration.

Usage: python3 -I h47_score.py <H47_dir> <extract_dir> <deg_results.csv> <out_dir>

Design is fixed in <H47_dir>/PRESPEC_H47.md. Scoring is the L1000-native weighted
connectivity score (WTCS) on level 5 landmark z-scores, normalized (NCS) within
cell line x series, summarized per compound (median within line, mean across
lines). Class recovery is an AUC / one-sided Mann-Whitney over compounds, BH
over the testable sets, with an expression-matched random gene set null.
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

SEED = 20260923
N_NULL = int(os.environ.get("H47_N_NULL", "1000"))  # override only for synthetic dry runs
N_BOOT = 2000
MIN_SIDE = 5  # PRESPEC addendum 1 (was 10)
MAX_SIDE = 150
TOP_T = 50
MIN_MEMBERS = 5
N_BINS = 10
LINES = ["HEPG2", "HUH7", "PHH"]
SERIES = ["GSE92742", "GSE70138"]

PPAR_TOKENS = {"ppar receptor agonist", "ppar receptor partial agonist"}
FXR_TOKENS = {"fxr agonist"}
THR_TOKENS = {"thyroid hormone receptor agonist"}
ACC_TOKENS = {"acetyl-coa carboxylase inhibitor"}

# Lefever et al. 2022, Metabolites 12:528, Table 1 (25 drugs), as printed.
N1_LEFEVER = [
    "eltanolone", "fenoprofen", "oxandrolone", "cefotaxime", "amorolfine",
    "dexamethasone", "proxyphylline", "sn-38", "sulfanitran", "tetracycline",
    "7-hydroxystaurosporine", "dopamine", "medrysone", "mestranol",
    "norethindrone", "troxerutin", "brequinar", "bromocriptine", "cebranopadol",
    "flucloxacillin", "granisetron", "hexestrol", "iohexol", "melphalan",
    "oxacillin",
]
# Khazaei-Poul et al. 2025, Gastroenterol Hepatol Bed Bench 18:230, Table 1
# plus lobelanidine (Results and Discussion), as printed.
N2_KHAZAEI = [
    "BRD-K59968781", "WYE-354", "BRD-K62176618", "BRD-K50617795", "AT-SUMO-1",
    "triciribine", "cinanserin", "fenoldopam", "BRD-K18710049", "apabetalone",
    "quinidine", "BRD-K12974237", "BRD-K53295358", "lobelanidine", "CGP-60474",
    "BRD-K15513102", "BRD-K26664453", "BRD-A59145032", "AZD-7762",
]
SETS = ["PPAR_agonist", "FXR_agonist", "N1_Lefever2022", "N2_KhazaeiPoul2025"]


def norm_name(s):
    return re.sub(r"[^0-9a-z]", "", str(s).casefold())


# ---------------------------------------------------------------- scoring core
class RankedRef:
    """Per-signature descending order of landmark z-scores, reused for every query."""

    def __init__(self, z):
        z = np.asarray(z, dtype=np.float64)
        self.order = np.argsort(-z, axis=1, kind="stable")
        self.absz = np.take_along_axis(np.abs(z), self.order, axis=1)
        self.n_sig, self.n_gene = z.shape

    def es(self, mask, chunk=4000):
        """Weighted (p = 1) KS enrichment score of gene set `mask` in every signature."""
        mask = np.asarray(mask, dtype=bool)
        k = int(mask.sum())
        out = np.empty(self.n_sig)
        for a in range(0, self.n_sig, chunk):
            b = min(a + chunk, self.n_sig)
            hits = mask[self.order[a:b]]
            w = self.absz[a:b] * hits
            nr = w.sum(axis=1)
            nr[nr == 0] = np.nan
            d = np.cumsum(w, axis=1) / nr[:, None] - np.cumsum(~hits, axis=1) / (self.n_gene - k)
            idx = np.nanargmax(np.abs(np.nan_to_num(d)), axis=1)
            out[a:b] = d[np.arange(b - a), idx]
        return np.nan_to_num(out)


def wtcs(es_up, es_dn):
    return np.where(np.sign(es_up) != np.sign(es_dn), (es_up - es_dn) / 2.0, 0.0)


def ncs(w, groups, n_groups):
    pos = w > 0
    neg = w < 0
    mu_pos = np.bincount(groups, weights=np.where(pos, w, 0.0), minlength=n_groups) / np.maximum(
        np.bincount(groups, weights=pos.astype(float), minlength=n_groups), 1)
    mu_neg = np.bincount(groups, weights=np.where(neg, -w, 0.0), minlength=n_groups) / np.maximum(
        np.bincount(groups, weights=neg.astype(float), minlength=n_groups), 1)
    scale = np.where(pos, mu_pos[groups], np.where(neg, mu_neg[groups], 1.0))
    return np.where(scale > 0, w / np.where(scale > 0, scale, 1.0), 0.0)


def compound_scores(r, comp, line, n_comp, n_line=len(LINES)):
    """Median r per compound x line, then mean over lines with data. Returns (pooled, per_line)."""
    key = comp * n_line + line
    med = pd.Series(r).groupby(key).median()
    per_line = np.full(n_comp * n_line, np.nan)
    per_line[med.index.to_numpy()] = med.to_numpy()
    per_line = per_line.reshape(n_comp, n_line)
    with np.errstate(invalid="ignore"):
        pooled = np.nanmean(per_line, axis=1)
    return pooled, per_line


def auc(scores, member):
    member = np.asarray(member, dtype=bool)
    n1 = member.sum()
    n0 = member.size - n1
    ranks = stats.rankdata(scores)
    return (ranks[member].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def mw_greater(scores, member):
    member = np.asarray(member, dtype=bool)
    return stats.mannwhitneyu(scores[member], scores[~member], alternative="greater",
                              method="asymptotic").pvalue


def boot_ci(scores, member, rng, n_boot=N_BOOT):
    member = np.asarray(member, dtype=bool)
    i1 = np.flatnonzero(member)
    i0 = np.flatnonzero(~member)
    lab = np.r_[np.ones(i1.size, bool), np.zeros(i0.size, bool)]
    vals = np.empty(n_boot)
    for b in range(n_boot):
        s = np.r_[scores[rng.choice(i1, i1.size)], scores[rng.choice(i0, i0.size)]]
        vals[b] = auc(s, lab)
    return np.percentile(vals, [2.5, 97.5])


def bh(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


# ---------------------------------------------------------------- inputs
def load_signatures(extract_dir):
    lm = pd.read_csv(os.path.join(extract_dir, "landmark_order.tsv"), sep="\t", dtype=str).pr_gene_id.tolist()
    zs, metas = [], []
    for s in SERIES:
        zs.append(np.load(os.path.join(extract_dir, f"z_landmark_{s}.npy")))
        m = pd.read_csv(os.path.join(extract_dir, f"sig_meta_{s}.tsv"), sep="\t", dtype=str)
        m["series"] = s
        metas.append(m)
    z = np.vstack(zs)
    meta = pd.concat(metas, ignore_index=True)
    assert z.shape[0] == len(meta) and z.shape[1] == len(lm)
    return z, meta, lm


def build_query(extract_dir, deg_path, lm):
    lmap = pd.read_csv(os.path.join(extract_dir, "landmark_map.tsv"), sep="\t", dtype=str).fillna("")
    lmap = lmap.set_index("pr_gene_id").loc[lm].reset_index()
    deg = pd.read_csv(deg_path)
    deg["ensembl"] = deg.gene.str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)
    deg = deg.set_index("ensembl")
    keep = (lmap.ensembl != "").to_numpy()
    rq = lmap[keep].copy()
    rq["col"] = np.flatnonzero(keep)
    rq = rq.join(deg[["logFC", "t", "padj", "AveExpr", "symbol"]], on="ensembl")
    assert rq.logFC.notna().all()
    rq["call"] = (rq.padj < 0.05) & (rq.logFC.abs() > 0.5)
    rq["direction"] = np.where(rq.call & (rq.logFC > 0), "up", np.where(rq.call & (rq.logFC < 0), "down", ""))
    for d in ["up", "down"]:
        idx = rq.index[rq.direction == d]
        if len(idx) > MAX_SIDE:
            drop = rq.loc[idx].logFC.abs().sort_values(ascending=False).index[MAX_SIDE:]
            rq.loc[drop, "direction"] = ""
    # Addendum 1 alternative query: top 50 by moderated t per side among padj < 0.05
    rq["dir_top50t"] = ""
    sig = rq[rq.padj < 0.05]
    rq.loc[sig[sig.t > 0].t.sort_values(ascending=False).index[:TOP_T], "dir_top50t"] = "up"
    rq.loc[sig[sig.t < 0].t.sort_values(ascending=True).index[:TOP_T], "dir_top50t"] = "down"
    rq["expr_bin"] = pd.qcut(rq.AveExpr.rank(method="first"), N_BINS, labels=False)
    return rq.reset_index(drop=True), int(keep.sum())


def compound_units(meta):
    pid = meta.pert_id.str.slice(0, 13)
    # GSE92742 name wins when the two series disagree for one pert_id
    name_by_pid = (meta.assign(pid=pid, nm=meta.pert_iname.str.casefold(),
                               pri=(meta.series != "GSE92742").astype(int))
                   .sort_values("pri").drop_duplicates("pid").set_index("pid").nm)
    unit = pid.map(name_by_pid)
    return unit, pid


def membership(units_df, dl_dir):
    """units_df: one row per compound unit with columns unit, pert_ids (set)."""
    drugs = pd.read_csv(os.path.join(dl_dir, "repurposing_drugs_20200324.txt"), sep="\t", comment="!", dtype=str).fillna("")
    samples = pd.read_csv(os.path.join(dl_dir, "repurposing_samples_20200324.txt"), sep="\t", comment="!", dtype=str).fillna("")
    cinfo = pd.read_csv(os.path.join(dl_dir, "compoundinfo_beta.txt"), sep="\t", dtype=str).fillna("")
    hub_tokens = {r.pert_iname: {t.strip().casefold() for t in r.moa.split("|") if t.strip()} for r in drugs.itertuples()}
    hub_name_cf = {k.casefold(): k for k in hub_tokens}
    pid2hub = samples.assign(pid=samples.broad_id.str.slice(0, 13)).groupby("pid").pert_iname.agg(set).to_dict()
    cinfo_tok = cinfo.assign(t=cinfo.moa.str.strip().str.casefold()).groupby("pert_id").t.agg(set).to_dict()
    name2pid = {}
    for r in cinfo.itertuples():
        for nm in (r.cmap_name, r.compound_aliases):
            if nm:
                name2pid.setdefault(norm_name(nm), set()).add(r.pert_id[:13])

    rows = []
    for r in units_df.itertuples():
        hub_names = set()
        for p in r.pert_ids:
            hub_names |= pid2hub.get(p, set())
        how = "broad_id"
        if not hub_names and r.unit in hub_name_cf:
            hub_names = {hub_name_cf[r.unit]}
            how = "name"
        if not hub_names:
            how = "none"
        toks = set().union(*[hub_tokens.get(h, set()) for h in hub_names]) if hub_names else set()
        ctoks = set().union(*[cinfo_tok.get(p, set()) for p in r.pert_ids]) if r.pert_ids else set()
        rows.append({
            "unit": r.unit, "hub_match": how, "hub_names": "|".join(sorted(hub_names)),
            "hub_moa": "|".join(sorted(toks)), "cinfo_moa": "|".join(sorted(t for t in ctoks if t)),
            "PPAR_agonist": bool(toks & PPAR_TOKENS), "FXR_agonist": bool(toks & FXR_TOKENS),
            "PPAR_agonist_cinfo": bool(ctoks & PPAR_TOKENS), "FXR_agonist_cinfo": bool(ctoks & FXR_TOKENS),
            "THR_agonist": bool(toks & THR_TOKENS), "ACC_inhibitor": bool(toks & ACC_TOKENS),
        })
    mem = pd.DataFrame(rows)

    # nominations
    unit_by_pid = {p: r.unit for r in units_df.itertuples() for p in r.pert_ids}
    unit_by_norm = {}
    for r in units_df.itertuples():
        unit_by_norm.setdefault(norm_name(r.unit), set()).add(r.unit)
    nom_rows = []
    for set_name, entries in [("N1_Lefever2022", N1_LEFEVER), ("N2_KhazaeiPoul2025", N2_KHAZAEI)]:
        hits_all = set()
        for e in entries:
            hits, how = set(), []
            if e.upper().startswith("BRD-"):
                if e.upper()[:13] in unit_by_pid:
                    hits.add(unit_by_pid[e.upper()[:13]])
                    how.append("pert_id")
            else:
                n = norm_name(e)
                if n in unit_by_norm:
                    hits |= unit_by_norm[n]
                    how.append("pert_iname")
                for p in name2pid.get(n, set()):
                    if p in unit_by_pid:
                        hits.add(unit_by_pid[p])
                        how.append("compoundinfo_alias")
            hits_all |= hits
            nom_rows.append({"set": set_name, "entry": e, "matched_units": "|".join(sorted(hits)),
                             "match_route": "|".join(sorted(set(how)))})
        mem[set_name] = mem.unit.isin(hits_all)
    return mem, pd.DataFrame(nom_rows)


# ---------------------------------------------------------------- pipeline
def run_pipeline(ref, up, dn, grp, n_grp, comp, line, n_comp, sig_keep=None):
    """dn=None gives the up-only score (ES_up, no WTCS sign rule)."""
    e_up = ref.es(up)
    if dn is None:
        e_dn = np.zeros_like(e_up)
        w = e_up
    else:
        e_dn = ref.es(dn)
        w = wtcs(e_up, e_dn)
    if sig_keep is None:
        sig_keep = np.ones(w.size, bool)
    n = np.zeros(w.size)
    n[sig_keep] = ncs(w[sig_keep], grp[sig_keep], n_grp)
    r = -n
    pooled, per_line = compound_scores(r[sig_keep], comp[sig_keep], line[sig_keep], n_comp)
    return dict(es_up=e_up, es_dn=e_dn, wtcs=w, ncs=n, r=r, pooled=pooled, per_line=per_line)


def set_tests(scores, mem, sets, rng, with_ci=True):
    ok = ~np.isnan(scores)
    out = []
    for s in sets:
        m = mem[s].to_numpy() & ok
        n1 = int(m.sum())
        row = {"set": s, "n_members": n1, "n_universe": int(ok.sum()), "testable": n1 >= MIN_MEMBERS}
        if n1 >= MIN_MEMBERS:
            sc = scores[ok]
            mm = m[ok]
            row["auc"] = auc(sc, mm)
            row["mw_p_one_sided"] = mw_greater(sc, mm)
            if with_ci:
                lo, hi = boot_ci(sc, mm, rng)
                row["auc_ci_lo"], row["auc_ci_hi"] = lo, hi
            pct = stats.rankdata(sc)[mm] / sc.size
            row["median_member_percentile"] = float(np.median(pct))
        out.append(row)
    return pd.DataFrame(out)


def matched_null(ref, up, dn, bins, grp, n_grp, comp, line, n_comp, mem, tests, label):
    """Expression-matched random gene sets with the query's per-bin up/down counts."""
    rng = np.random.default_rng(SEED)
    bin_genes = [np.flatnonzero(bins == b) for b in range(N_BINS)]
    n_up_b = [int(up[g].sum()) for g in bin_genes]
    n_dn_b = [int(dn[g].sum()) if dn is not None else 0 for g in bin_genes]
    testable = [s for s in tests.loc[tests.testable, "set"]]
    rows = []
    for d in range(N_NULL):
        u = np.zeros(up.size, bool)
        v = np.zeros(up.size, bool)
        for g, nu, nd in zip(bin_genes, n_up_b, n_dn_b):
            pick = rng.choice(g, nu + nd, replace=False)
            u[pick[:nu]] = True
            v[pick[nu:]] = True
        res = run_pipeline(ref, u, v if dn is not None else None, grp, n_grp, comp, line, n_comp)
        for s in testable:
            m = mem[s].to_numpy()
            rows.append({"query": label, "draw": d, "set": s, "auc": auc(res["pooled"], m),
                         "mw_p_one_sided": mw_greater(res["pooled"], m)})
        if d % 100 == 0:
            print(f"[{label}] null draw {d}", flush=True)
    null = pd.DataFrame(rows, columns=["query", "draw", "set", "auc", "mw_p_one_sided"])
    tests = tests.copy()
    for i, row in tests.iterrows():
        if not row.testable:
            continue
        nd = null[null.set == row.set]
        tests.loc[i, "p_null"] = (1 + (nd.auc >= row.auc).sum()) / (1 + len(nd))
        tests.loc[i, "null_auc_median"] = nd.auc.median()
        tests.loc[i, "null_auc_q025"] = nd.auc.quantile(0.025)
        tests.loc[i, "null_auc_q975"] = nd.auc.quantile(0.975)
        tests.loc[i, "null_frac_mw_p_lt_0.05"] = (nd.mw_p_one_sided < 0.05).mean()
    return tests, null


def main(h47_dir, extract_dir, deg_path, out_dir):
    os.makedirs(out_dir, exist_ok=False)
    dl_dir = os.path.join(h47_dir, "downloads_untrusted")
    summary = {"seed": SEED, "n_null": N_NULL, "n_boot": N_BOOT, "min_side": MIN_SIDE}

    z, meta, lm = load_signatures(extract_dir)
    rq, n_mapped = build_query(extract_dir, deg_path, lm)
    summary["landmarks_mapped"] = n_mapped
    summary["landmark_map_flag_lt900"] = n_mapped < 900
    rq.to_csv(os.path.join(out_dir, "ranked_universe_and_query.tsv"), sep="\t", index=False)
    up = (rq.direction == "up").to_numpy()
    dn = (rq.direction == "down").to_numpy()
    up_t = (rq.dir_top50t == "up").to_numpy()
    dn_t = (rq.dir_top50t == "down").to_numpy()
    summary.update(query_up=int(up.sum()), query_down=int(dn.sum()),
                   landmark_padj05_up=int(((rq.padj < 0.05) & (rq.logFC > 0)).sum()),
                   landmark_padj05_down=int(((rq.padj < 0.05) & (rq.logFC < 0)).sum()),
                   top50t_up=int(up_t.sum()), top50t_down=int(dn_t.sum()),
                   top50t_overlap_with_canonical=int(((up_t & up) | (dn_t & dn)).sum()))
    if up.sum() < MIN_SIDE or dn.sum() < MIN_SIDE:
        summary["stopped"] = f"S3: query up {int(up.sum())}, down {int(dn.sum())} (< {MIN_SIDE})"
        json.dump(summary, open(os.path.join(out_dir, "summary.json"), "w"), indent=2)
        sys.exit("STOP " + summary["stopped"])

    ref = RankedRef(z[:, rq.col.to_numpy()])
    bins = rq.expr_bin.to_numpy()

    unit, pid = compound_units(meta)
    meta = meta.assign(unit=unit.to_numpy(), pid=pid.to_numpy())
    units = meta.groupby("unit").pid.agg(lambda x: set(x)).reset_index().rename(columns={"pid": "pert_ids"})
    mem, nom = membership(units, dl_dir)
    comp_codes = {u: i for i, u in enumerate(mem.unit)}
    comp = meta.unit.map(comp_codes).to_numpy()
    line = meta.cell_id.map({l: i for i, l in enumerate(LINES)}).to_numpy()
    assert line.dtype.kind == "i"
    grp_key = meta.cell_id + "|" + meta.series
    grp_codes = {g: i for i, g in enumerate(sorted(grp_key.unique()))}
    grp = grp_key.map(grp_codes).to_numpy()
    n_comp, n_grp = len(mem), len(grp_codes)
    args = (grp, n_grp, comp, line, n_comp)
    nom.to_csv(os.path.join(out_dir, "nomination_matching.tsv"), sep="\t", index=False)
    summary["n_signatures"] = int(len(meta))
    summary["signatures_per_line_series"] = grp_key.value_counts().to_dict()
    summary["n_compounds"] = n_comp
    summary["compounds_per_line"] = {l: int(meta.loc[meta.cell_id == l, "unit"].nunique()) for l in LINES}
    summary["members"] = {s: int(mem[s].sum()) for s in SETS + ["PPAR_agonist_cinfo", "FXR_agonist_cinfo", "THR_agonist", "ACC_inhibitor"]}
    summary["hub_match_routes"] = mem.hub_match.value_counts().to_dict()

    # ---- observed, three queries
    obs = run_pipeline(ref, up, dn, *args)
    obs_t = run_pipeline(ref, up_t, dn_t, *args)
    obs_u = run_pipeline(ref, up, None, *args)
    sig_out = meta[["sig_id", "pert_id", "pert_iname", "unit", "cell_id", "series", "pert_idose", "pert_itime",
                    "distil_cc_q75", "distil_nsample"]].copy()
    for k in ["es_up", "es_dn", "wtcs", "ncs", "r"]:
        sig_out[k] = obs[k]
    sig_out["r_top50t"] = obs_t["r"]
    sig_out["r_up_only"] = obs_u["r"]
    sig_out.to_csv(os.path.join(out_dir, "signature_scores.tsv.gz"), sep="\t", index=False)

    comp_out = mem.copy()
    comp_out["pooled_reversal"] = obs["pooled"]
    for i, l in enumerate(LINES):
        comp_out[f"reversal_{l}"] = obs["per_line"][:, i]
        comp_out[f"n_sig_{l}"] = np.bincount(comp[line == i], minlength=n_comp)
    comp_out["rank"] = stats.rankdata(-obs["pooled"], method="min").astype(int)
    comp_out["percentile"] = stats.rankdata(obs["pooled"]) / n_comp
    comp_out["pooled_reversal_top50t"] = obs_t["pooled"]
    comp_out["pooled_reversal_up_only"] = obs_u["pooled"]
    comp_out.sort_values("rank").to_csv(os.path.join(out_dir, "compound_scores.tsv"), sep="\t", index=False)

    rng_boot = np.random.default_rng(SEED)
    prim = set_tests(obs["pooled"], mem, SETS, rng_boot)
    t = prim.testable.to_numpy()
    prim["bh_q"] = np.nan
    if t.any():
        prim.loc[t, "bh_q"] = bh(prim.loc[t, "mw_p_one_sided"].to_numpy())
    prim, null_p = matched_null(ref, up, dn, bins, *args, mem, prim, "primary_canonical")
    if "p_null" not in prim:
        prim["p_null"] = np.nan
    prim["recovered"] = prim.testable & (prim.bh_q < 0.05) & (prim.p_null < 0.05)
    prim.insert(0, "query", "primary_canonical")
    prim.to_csv(os.path.join(out_dir, "class_tests_primary.tsv"), sep="\t", index=False)

    # ---- alternative queries (addendum 1; descriptive, own null)
    alt, nulls = [], [null_p]
    for label, res, u, v in [("top50_t", obs_t, up_t, dn_t), ("up_only", obs_u, up, None)]:
        tt = set_tests(res["pooled"], mem, SETS, np.random.default_rng(SEED))
        tt, nl = matched_null(ref, u, v, bins, *args, mem, tt, label)
        tt.insert(0, "query", label)
        alt.append(tt)
        nulls.append(nl)
    pd.concat(alt).to_csv(os.path.join(out_dir, "class_tests_alt_queries.tsv"), sep="\t", index=False)
    pd.concat(nulls).to_csv(os.path.join(out_dir, "null_draws.tsv.gz"), sep="\t", index=False)

    # ---- descriptive: per line, and HEPG2 per series (primary query)
    desc = []
    for i, l in enumerate(LINES):
        d = set_tests(obs["per_line"][:, i], mem, SETS, rng_boot)
        d.insert(0, "scope", l)
        desc.append(d)
    for s in SERIES:
        keep = ((meta.cell_id == "HEPG2") & (meta.series == s)).to_numpy()
        pooled, _ = compound_scores(obs["r"][keep], comp[keep], line[keep], n_comp)
        d = set_tests(pooled, mem, SETS, rng_boot)
        d.insert(0, "scope", f"HEPG2_{s}")
        desc.append(d)
    pd.concat(desc).to_csv(os.path.join(out_dir, "class_tests_per_line.tsv"), sep="\t", index=False)

    # ---- sensitivity (primary query, descriptive)
    sens = []
    cc = pd.to_numeric(meta.distil_cc_q75, errors="coerce").to_numpy()
    for name, keep in [("reproducible_cc_q75_ge_0.2", cc >= 0.2),
                       ("GSE92742_only", (meta.series == "GSE92742").to_numpy())]:
        res = run_pipeline(ref, up, dn, *args, sig_keep=keep)
        d = set_tests(res["pooled"], mem, SETS, rng_boot)
        d.insert(0, "analysis", name)
        d["n_signatures"] = int(keep.sum())
        sens.append(d)
    mem2 = mem.assign(PPAR_agonist=mem.PPAR_agonist_cinfo, FXR_agonist=mem.FXR_agonist_cinfo)
    d = set_tests(obs["pooled"], mem2, ["PPAR_agonist", "FXR_agonist"], rng_boot)
    d.insert(0, "analysis", "annotation_LINCS2020_compoundinfo")
    d["n_signatures"] = int(len(meta))
    sens.append(d)
    pd.concat(sens).to_csv(os.path.join(out_dir, "class_tests_sensitivity.tsv"), sep="\t", index=False)

    # ---- member ranks (primary query)
    keep_cols = (["unit", "hub_moa", "pooled_reversal", "rank", "percentile"] + [f"reversal_{l}" for l in LINES]
                 + [f"n_sig_{l}" for l in LINES] + ["pooled_reversal_top50t", "pooled_reversal_up_only"])
    mr = []
    for s in SETS + ["PPAR_agonist_cinfo", "FXR_agonist_cinfo", "THR_agonist", "ACC_inhibitor"]:
        x = comp_out.loc[comp_out[s], keep_cols].copy()
        x.insert(0, "set", s)
        mr.append(x)
    pd.concat(mr).sort_values(["set", "rank"]).to_csv(os.path.join(out_dir, "member_ranks.tsv"), sep="\t", index=False)

    summary["primary"] = json.loads(prim.to_json(orient="records"))
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(prim.to_string())


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
