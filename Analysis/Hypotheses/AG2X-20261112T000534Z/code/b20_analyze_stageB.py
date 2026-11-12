#!/usr/bin/env python3
"""Stage B analysis exactly as fixed in PRESPEC_AG2X_stageB.md (sections 4-5).

No model, calibration or prior is fitted on AlphaGenome outputs: only percentiles against the scored
background, permutation nulls and counts. B3 exact Mantel-Haenszel tests run in b21_b3_mh.R on the table
written here. Unit of inference: 1-Mb block chr~floor(pos_hg38/1e6). Seed 20261111. Nothing is adopted.
"""
from __future__ import annotations

import gzip
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

B = L.PKG / "stageB"
RAW, IN, TAB = B / "raw", B / "inputs", B / "tables"
FM = L.REPO / "GWAS/finemapping/results"
CURRIN = L.REPO / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1"
D1 = FM / "alphagenome_program/d1-lineage-allelic-20260914T194213Z/tables/lineage_site_summary.tsv"
MPRA = FM / "seqfunc/mpra_benchmark/mpra_substrate_truth.tsv"
NOM = FM / "seqfunc/nominated_cs_members.tsv"
N_PERM = 10_000
SEED = L.SEED
CONDS = ["bulk_liver", "hepatocyte", "monocyte_ctrl", "neuron_ctrl", "stellate", "kupffer", "liver_endothelial",
         "cholangiocyte", "hepg2"]
BENCH = {"BL_ATAC": "bulk_liver|ATAC", "BL_DNASE": "bulk_liver|DNASE", "HEP_SCATAC": "hepatocyte|SC_ATAC"}
CTRL_SCATAC = ["monocyte_ctrl|SC_ATAC", "neuron_ctrl|SC_ATAC"]
NOTES = {"candidate_not_adopted": True, "package_version": L.package_version(), "seed": SEED}


def block(uid):
    c, p = uid.split(":")[:2]
    return f"{c}~{int(p) // 1_000_000}"


def akey(k):
    c, p, a, b = k.split(":")
    return ":".join([c, p, min(a, b), max(a, b)])


def flip(uid):
    c, p, r, a = uid.split(":")
    return f"{c}:{p}:{a}:{r}"


# ------------------------------------------------------------------ scores
def load_scores():
    cols = ["variant_id", "output_type", "track_name", "imputation_mode", "raw_score"]
    vals, tracks = [], []
    for d in sorted(RAW.glob("t?_*")):
        cond = d.name.split("_", 1)[1]
        files = sorted(d.glob("chunk_*.parquet"))
        if not files:
            continue
        x = pd.concat([pd.read_parquet(f, columns=cols) for f in files], ignore_index=True)
        x = x[x.output_type.isin(["ATAC", "DNASE", "SC_ATAC"])]
        t = x.groupby(["output_type", "imputation_mode", "track_name"]).size().rename("rows").reset_index()
        t.insert(0, "condition", cond)
        t.insert(0, "tier_dir", d.name)
        tracks.append(t)
        x = x[x.imputation_mode == "observed"]
        x["uid"] = x.variant_id.str.replace(">", ":", regex=False)
        g = x.groupby(["uid", "output_type"]).raw_score.mean().reset_index()
        g["readout"] = cond + "|" + g.output_type
        vals.append(g[["uid", "readout", "raw_score"]])
        L.log(f"{d.name}: {x.uid.nunique()} variants with observed rows")
    tracks = pd.concat(tracks, ignore_index=True)
    tracks.to_csv(TAB / "stageB_observed_tracks_used.tsv", sep="\t", index=False)
    s = pd.concat(vals, ignore_index=True).drop_duplicates(["uid", "readout"])
    return s.pivot(index="uid", columns="readout", values="raw_score")


def percentiles(raw, bg_uids):
    """Share of the pooled background |readout| strictly below each variant's |readout| (per readout)."""
    pct = pd.DataFrame(index=raw.index, columns=raw.columns, dtype=float)
    thr = {}
    for r in raw.columns:
        b = np.sort(np.abs(raw.loc[raw.index.intersection(bg_uids), r].dropna().to_numpy()))
        thr[r] = {"n_background": int(len(b)), "abs_p99": float(np.quantile(b, 0.99)) if len(b) else None}
        x = raw[r].abs()
        ok = x.notna()
        pct.loc[ok, r] = np.searchsorted(b, x[ok].to_numpy(), side="left") / len(b)
    return pct, thr


# ------------------------------------------------------------------ B1
def b1(sets, pct, rng_seed):
    """sets: rows (set_id, uid, pip, block). Statistic T = mean over blocks of mean over sets of PIP-weighted pct."""
    out = []
    readouts = {"B1a_HEP_SCATAC": pct[BENCH["HEP_SCATAC"]],
                "B1b_HEP_minus_controls": pct[BENCH["HEP_SCATAC"]] - pct[CTRL_SCATAC].mean(axis=1, skipna=False),
                "BL_ATAC": pct[BENCH["BL_ATAC"]], "BL_DNASE": pct[BENCH["BL_DNASE"]]}
    for r in pct.columns:
        if r.split("|")[0] in ("monocyte_ctrl", "neuron_ctrl"):
            readouts[f"control_{r}"] = pct[r]
    for name, val in readouts.items():
        d = sets.assign(v=sets.uid.map(val)).dropna(subset=["v"])
        d = d[d.groupby("set_id").uid.transform("size") >= 2].sort_values(["set_id", "uid"]).reset_index(drop=True)
        if d.empty:
            out.append({"test": name, "n_sets": 0})
            continue
        si = pd.factorize(d.set_id)[0]
        set_block = d.groupby(si).block.first().to_numpy()
        bi, blocks = pd.factorize(set_block)
        pip, v = d.pip.to_numpy(), d.v.to_numpy()
        spip = np.bincount(si, pip)
        nb_sets = np.bincount(bi)

        def stat(p):
            w = np.bincount(si, p * v) / spip
            return float(np.mean(np.bincount(bi, w) / nb_sets))

        t_obs = stat(pip)
        rng = np.random.default_rng(rng_seed)  # same draws for every readout
        ge = 0
        null = np.empty(N_PERM)
        for i in range(N_PERM):
            order = np.lexsort((rng.random(len(si)), si))
            null[i] = stat(pip[order])
            ge += null[i] >= t_obs
        unweighted = float(np.mean(np.bincount(bi, np.bincount(si, v) / np.bincount(si)) / nb_sets))
        out.append({"test": name, "n_sets": int(si.max() + 1), "n_members": int(len(d)), "n_blocks": int(len(blocks)),
                    "T_obs": t_obs, "T_unweighted": unweighted, "null_mean": float(null.mean()),
                    "null_q025": float(np.quantile(null, 0.025)), "null_q975": float(np.quantile(null, 0.975)),
                    "p_one_sided": (1 + ge) / (N_PERM + 1)})
    res = pd.DataFrame(out)
    fam = res.test.isin(["B1a_HEP_SCATAC", "B1b_HEP_minus_controls"])
    p = res.loc[fam, "p_one_sided"].to_numpy()
    order = np.argsort(p)
    holm = np.empty_like(p)
    running = 0.0
    for rank, j in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[j]))
        holm[j] = running
    res["p_holm"] = np.nan
    res.loc[fam, "p_holm"] = holm
    res["decision"] = np.where(fam, np.where(res.p_holm < 0.05, "pass", "fail"), "reported, no test")
    return res


# ------------------------------------------------------------------ B2
def b2(loci_members, loci_bg, supported):
    rng = np.random.default_rng(SEED)
    obs, pools = 0, []
    for lid, mem in loci_members.items():
        m = [supported.get(u, np.nan) for u in mem]
        b = [supported.get(u, np.nan) for u in loci_bg.get(lid, [])]
        m = [x for x in m if x == x]
        b = [x for x in b if x == x]
        if not m:
            continue
        obs += int(sum(m) == 1)
        pools.append((np.array(m + b, dtype=bool), len(m)))
    null = np.zeros(N_PERM, dtype=int)
    for i in range(N_PERM):
        c = 0
        for pool, k in pools:
            c += int(pool[rng.permutation(len(pool))[:k]].sum() == 1)
        null[i] = c
    return {"n_loci": len(pools), "observed_exactly_one": obs, "null_mean": float(null.mean()),
            "null_q025": float(np.quantile(null, 0.025)), "null_q975": float(np.quantile(null, 0.975)),
            "p_upper": float((1 + (null >= obs).sum()) / (N_PERM + 1))}


# ------------------------------------------------------------------ B3 outcomes
def caqtl_effects(uids):
    """ALT-allele caQTL beta per scored uid (largest |beta| if several peaks). Proxy rule per prespec."""
    lead = pd.read_csv(CURRIN / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
                       sep="\t", usecols=["lead_variant_ID", "peak_ID", "beta", "EA", "NEA"])
    lc = lead.lead_variant_ID.str.split(":", expand=True)
    lead["beta_alt"] = np.where(lead.EA == lc[3], lead.beta, np.where(lead.EA == lc[2], -lead.beta, np.nan))
    prox = pd.read_csv(CURRIN / "caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz", sep="\t",
                       usecols=["proxy_ID", "lead_ID", "caPeak", "proxy_caQTL_beta"])
    chk = prox[prox.proxy_ID == prox.lead_ID].merge(lead, left_on=["lead_ID", "caPeak"],
                                                    right_on=["lead_variant_ID", "peak_ID"])
    chk = chk.dropna(subset=["beta_alt"])
    agree = float((np.sign(chk.proxy_caQTL_beta) == np.sign(chk.beta_alt)).mean()) if len(chk) else float("nan")
    use_proxies = bool(agree >= 0.99)
    NOTES["B3a_proxy_orientation_check"] = {"lead_rows_compared": int(len(chk)), "sign_agreement": agree,
                                            "proxies_used": use_proxies}
    eff = lead[["lead_variant_ID", "beta_alt"]].rename(columns={"lead_variant_ID": "id"}).dropna()
    if use_proxies:
        eff = pd.concat([eff, prox.rename(columns={"proxy_ID": "id", "proxy_caQTL_beta": "beta_alt"})[["id", "beta_alt"]]])
    eff = eff.assign(a=eff.beta_alt.abs()).sort_values("a", ascending=False).drop_duplicates("id")
    m = dict(zip(eff.id, eff.beta_alt))
    out = {}
    for u in uids:
        if u in m:
            out[u] = m[u]
        elif flip(u) in m:
            out[u] = -m[flip(u)]
    NOTES["B3a_members_with_caqtl"] = len(out)
    return out


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    v = pd.read_csv(IN / "variants_scored_final.tsv", sep="\t")
    imap_path = IN / "indel_uid_map.tsv"
    imap = dict(pd.read_csv(imap_path, sep="\t").values) if imap_path.exists() else {}
    fix = lambda u: imap.get(u, u)  # noqa: E731
    raw = load_scores()
    bg = set(v.loc[v.roles.str.contains("background"), "uid"])
    pct, thr = percentiles(raw, bg)
    NOTES["thresholds"] = thr
    NOTES["variants_with_scores"] = int(raw.shape[0])

    # ---------------- per-variant fields
    mv = pd.read_csv(IN / "member_variants.tsv", sep="\t", dtype={"coding_pip05": str})
    mv["uid"] = mv.uid.map(fix)
    coding = dict(zip(mv.uid, mv.coding_pip05.str.upper() == "TRUE"))
    pv = v.set_index("uid")[["tier", "roles", "loci"]].join(raw.add_prefix("raw:")).join(pct.add_prefix("pct:"))
    sup = (pct >= 0.99)
    sup = sup.where(pct.notna())
    for k, r in BENCH.items():
        pv[f"supported_{k}"] = sup[r]
    is_member = pv.roles.str.contains(r"(^|;)member($|;)")
    bench_sup = sup[list(BENCH.values())].reindex(pv.index)
    pv["sequence_supported"] = is_member & bench_sup.fillna(False).astype(bool).any(axis=1)
    cell = sup.reindex(pv.index)
    pv["supported_conditions"] = cell.apply(
        lambda row: ";".join(sorted({c.split("|")[0] + ("(benchmarked)" if c in BENCH.values() else "(not benchmarked)")
                                     for c, x in row.items() if x == 1.0})), axis=1)

    def acc_sign(u):
        s = {np.sign(raw.at[u, r]) for k, r in BENCH.items() if sup.at[u, r] == 1.0} if u in raw.index else set()
        return "" if not s else ("gain" if s == {1.0} else "loss" if s == {-1.0} else "mixed")
    pv["accessibility_sign_alt_minus_ref"] = [acc_sign(u) for u in pv.index]
    pv["coding_pip05"] = pv.index.map(lambda u: coding.get(u, False))
    pv["candidate_not_adopted"] = True
    pv["package_version"] = NOTES["package_version"]
    pv.reset_index().to_csv(TAB / "stageB_variant_fields.tsv.gz", sep="\t", index=False, compression="gzip")

    # ---------------- per-locus fields (loci sorted by ID; no rank across loci)
    loci = pd.read_csv(IN / "loci.tsv", sep="\t")
    rows = []
    memb = pv[is_member]
    for _, lc in loci.sort_values("locus_id").iterrows():
        lid = lc.locus_id
        here = memb[memb.loci.str.split(";").map(lambda x: lid in x)]
        if lc.scored_as != "credible_set_members":
            field = "not applicable (scored by lead and proxies)" if lc.scored_as != "not_scored" else "not scored"
            rows.append({"locus_id": lid, "tier": lc.tier, "scored_as": lc.scored_as, "sequence_supported_variant": field})
            continue
        ss = here[here.sequence_supported]
        n = len(ss)
        rows.append({"locus_id": lid, "tier": lc.tier, "genes": lc.get("genes", ""), "scored_as": lc.scored_as,
                     "n_members_scored": int(here[f"raw:{BENCH['BL_ATAC']}"].notna().sum()) if len(here) else 0,
                     "sequence_supported_variant": "none (not evidence against the locus)" if n == 0 else
                     ("one" if n == 1 else "several"),
                     "supported_uids": ";".join(ss.index), "supported_conditions": "|".join(ss.supported_conditions),
                     "accessibility_sign": ";".join(ss.accessibility_sign_alt_minus_ref),
                     "coding_members_pip05": ";".join(here.index[here.coding_pip05])})
    lt = pd.DataFrame(rows)
    lt["candidate_not_adopted"] = True
    lt.to_csv(TAB / "stageB_locus_fields.tsv", sep="\t", index=False)

    # ---------------- B1
    with gzip.open(IN / "credible_set_member_rows.tsv.gz", "rt") as fh:
        cs = pd.read_csv(fh, sep="\t", usecols=["locus_id", "set_id", "k", "method", "pip", "uid", "tier"])
    cs["uid"] = cs.uid.map(lambda u: fix(u) if isinstance(u, str) else u)
    cs = cs[cs.uid.isin(raw.index)]
    su = cs[cs.method == "susie"].drop_duplicates(["set_id", "uid"])
    top = su.sort_values(["set_id", "pip"], ascending=[True, False]).drop_duplicates("set_id")
    su = su.merge(top[["set_id", "uid"]].rename(columns={"uid": "top_uid"}), on="set_id")
    su["block"] = su.top_uid.map(block)
    b1res = b1(su[["set_id", "uid", "pip", "block"]], pct, SEED)
    b1res["candidate_not_adopted"] = True
    b1res.to_csv(TAB / "stageB_B1.tsv", sep="\t", index=False)

    # ---------------- B2
    cs_loci = set(loci.loc[loci.scored_as == "credible_set_members", "locus_id"])
    loci_members = {lid: list(set(g.uid)) for lid, g in cs[cs.locus_id.isin(cs_loci)].groupby("locus_id")}
    bgt = pd.read_csv(IN / "background.tsv", sep="\t")
    loci_bg = {lid: list(set(g.uid)) for lid, g in bgt.groupby("locus_id")}
    seqsup = {u: float(bool(x)) for u, x in bench_sup.fillna(False).astype(bool).any(axis=1).items()
              if u in raw.index}
    r2 = b2(loci_members, loci_bg, seqsup)
    nom = pd.read_csv(NOM, sep="\t", usecols=["variant_id"])
    nomk = set(nom.variant_id.map(akey))
    k_of = dict(zip(mv.uid, mv.k.map(akey)))
    ss_uids = pv.index[pv.sequence_supported]
    r2["sequence_supported_members"] = int(len(ss_uids))
    r2["of_which_in_nominated_cs_members"] = int(sum(k_of.get(u) in nomk for u in ss_uids))
    r2["nominated_members_scored"] = int(sum(k in nomk for k in k_of.values()))
    r2["candidate_not_adopted"] = True
    json.dump(r2, open(TAB / "stageB_B2.json", "w"), indent=1)

    # ---------------- B3 input (sets with >= 1 supported and >= 1 unsupported member, per arm)
    mem_cs = cs[["locus_id", "set_id", "uid"]].drop_duplicates()
    ca = caqtl_effects(mem_cs.uid.unique())
    d1 = pd.read_csv(D1, sep="\t")
    d1 = d1[(d1.scope == "POOLED_30_DONORS") & (d1.lineage == "Hepatocyte")]
    ai = {u: float(abs(m / s) >= 1.96) for u, m, s in zip(d1.uid, d1.mean_log2_alt_over_ref, d1.se_log2) if s > 0}
    ai.update({flip(u): x for u, x in list(ai.items()) if flip(u) not in ai})
    mp = pd.read_csv(MPRA, sep="\t", dtype=str)
    mp = mp[mp.context.isin(["HepG2_ctrl", "HepG2_PAOA"]) & (mp.mpra_barcode_mapped.str.upper() == "TRUE")]
    mp["uid"] = mp.chr + ":" + mp.pos_hg38 + ":" + mp.ref_hg38 + ":" + mp.alt_hg38
    mpd = (mp.assign(d=mp.dav.astype(float)).groupby("uid").d.max()).to_dict()
    mpd.update({flip(u): x for u, x in list(mpd.items()) if flip(u) not in mpd})
    arms = {"B3a_caqtl_matching_sign": (sup[BENCH["BL_ATAC"]],
                                        lambda u: float(u in ca and np.sign(ca[u]) == np.sign(raw.at[u, BENCH["BL_ATAC"]]))),
            "B3b_hepatocyte_allelic_imbalance": (sup[BENCH["HEP_SCATAC"]], lambda u: ai.get(u, np.nan)),
            "B3c_mpra_allelic_activity": (pd.Series(seqsup), lambda u: mpd.get(u, np.nan))}
    b3rows = []
    for arm, (s, outcome) in arms.items():
        d = mem_cs.assign(supported=mem_cs.uid.map(s)).dropna(subset=["supported"])
        ok = d.groupby("set_id").supported.agg(lambda x: x.min() == 0 and x.max() == 1)
        d = d[d.set_id.isin(ok[ok].index)].drop_duplicates(["locus_id", "uid"])
        d["outcome"] = d.uid.map(outcome)
        d = d.dropna(subset=["outcome"])
        d["arm"] = arm
        b3rows.append(d[["arm", "locus_id", "uid", "supported", "outcome"]])
    pd.concat(b3rows).to_csv(TAB / "stageB_B3_input.tsv", sep="\t", index=False)
    a1 = pd.read_csv(L.REPO / "GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z/tables/"
                     "endpoint1_matched_tierA4.tsv.gz", sep="\t", usecols=["chr_x", "pos_hg38", "ref", "alt"], dtype=str)
    a1c = a1.chr_x.where(a1.chr_x.str.startswith("chr"), "chr" + a1.chr_x)
    a1u = set(a1c + ":" + a1.pos_hg38 + ":" + a1.ref + ":" + a1.alt)
    a1u |= {flip(u) for u in a1u}
    NOTES["B3a_members_in_stageA_A1"] = int(len(set(ca) & a1u))
    NOTES["B3b_members_in_stageA_A2_sites"] = int(len(set(mem_cs.uid) & set(d1.uid)))

    # ---------------- B4 (descriptive)
    dl = loci[(loci.source == "DIG1") & (loci.scored_as == "credible_set_members")]
    hep_sup = set(pv.index[is_member & (pv.supported_HEP_SCATAC == 1.0)])
    b4 = []
    for eff, g in dl.groupby("disease_effect"):
        n_with = sum(any(u in hep_sup for u in loci_members.get(lid, [])) for lid in g.locus_id)
        b4.append({"disease_effect": bool(eff), "n_loci": int(len(g)), "n_with_HEP_SCATAC_supported_member": int(n_with),
                   "share": n_with / len(g) if len(g) else np.nan, "candidate_not_adopted": True})
    pd.DataFrame(b4).to_csv(TAB / "stageB_B4.tsv", sep="\t", index=False)

    # ---------------- controls
    ctl = {}
    s1 = "chr1:109274968:G:T"
    if s1 in raw.index:
        bench = {k: {"raw": float(raw.at[s1, r]), "pct": float(pct.at[s1, r])} for k, r in BENCH.items()
                 if pd.notna(raw.at[s1, r])}
        cpct = [float(pct.at[s1, r]) for r in pct.columns if r.split("|")[0] in ("monocyte_ctrl", "neuron_ctrl")
                and pd.notna(pct.at[s1, r])]
        met = any(b["raw"] > 0 and b["pct"] > max(cpct) for b in bench.values()) if cpct else None
        ctl["SORT1_rs12740374"] = {"benchmarked": bench, "control_percentiles": cpct, "expectation_met": met}
    sp = sorted((RAW / "control_hsd17b13_splice").glob("chunk_*.parquet"))
    if sp:
        x = pd.concat([pd.read_parquet(f) for f in sp], ignore_index=True)
        ctl["HSD17B13_rs72613567_splice"] = (x.groupby("output_type").raw_score
                                             .agg(["size", lambda s: float(s.abs().max())]).reset_index()
                                             .rename(columns={"size": "rows", "<lambda_0>": "max_abs_raw"})
                                             .to_dict(orient="records"))
    ctl["coding_controls"] = {u: bool(coding.get(u, False)) for u in
                              ["chr22:43928847:C:G", "chr19:19268740:C:T", "chr2:27508073:T:C"]}
    json.dump(ctl, open(TAB / "stageB_controls.json", "w"), indent=1, default=str)

    # ---------------- repeat check
    cols = ["variant_id", "output_type", "track_name", "raw_score"]
    rp = sorted((RAW / "repeat_bulk_liver").glob("chunk_*.parquet"))
    if rp:
        a = pd.concat([pd.read_parquet(f, columns=cols) for f in rp], ignore_index=True)
        b = pd.concat([pd.read_parquet(f, columns=cols) for f in sorted((RAW / "t1_bulk_liver").glob("chunk_*.parquet"))],
                      ignore_index=True)
        m = a.merge(b, on=["variant_id", "output_type", "track_name"], suffixes=("_rep", "_orig"))
        diff = (m.raw_score_rep - m.raw_score_orig).abs()
        NOTES["repeat"] = {"rows_compared": int(len(m)), "variants": int(m.variant_id.nunique()),
                           "share_identical": float((diff == 0).mean()), "max_abs_diff": float(diff.max())}
    json.dump(NOTES, open(TAB / "stageB_notes.json", "w"), indent=1, default=str)
    L.log(json.dumps({k: NOTES[k] for k in NOTES if k != "thresholds"}, default=str))


if __name__ == "__main__":
    main()
