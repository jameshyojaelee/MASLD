#!/usr/bin/env python3
"""Stage B2 post-review 2026-11-18 analysis on cached parquets (zero API calls).

Writes exclusively to out_20261118/stageB2/ without touching stageB2/tables*/.
Applies the same N06 fixes as b20_analyze_stageB_20261118.py:
  1. Calibrated 3-track union threshold q99_union, empirical FDR, BH q, max_pip >= 0.10 / >= 0.50.
  2. Locus-grouped within_set_perm in B1 (unique uid per physical locus permuted once).
  3. Returns np.nan (NOT 0.0) in B3a when uid is unmeasured in Currin.
"""
from __future__ import annotations

import gzip
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402
import b20_analyze_stageB_20261118 as S  # noqa: E402

B2 = L.PKG / "stageB2"
RAW, IN = B2 / "raw", B2 / "inputs"
RAWB = L.PKG / "stageB/raw"
OUT = L.PKG / "out_20261118/stageB2"
BENCH = S.BENCH
NOTES = {
    "candidate_not_adopted": True,
    "package_version": L.package_version(),
    "seed": L.SEED,
    "model_version": "ModelVersion.ALL_FOLDS",
    "training_exposure_status": "ALL_FOLDS_TRAINING_EXPOSED_ZERO_GENOMIC_HOLDOUT",
}


def flip(u: str) -> str:
    c, p, r, a = u.split(":")
    return f"{c}:{p}:{a}:{r}"


def load_dirs(dirs, keep_uids=None):
    cols = ["variant_id", "output_type", "track_name", "imputation_mode", "raw_score"]
    vals, tracks = [], []
    for d, cond in dirs:
        files = sorted(d.glob("chunk_*.parquet"))
        if not files:
            continue
        x = pd.concat([pd.read_parquet(f, columns=cols) for f in files], ignore_index=True)
        x = x[x.output_type.isin(["ATAC", "DNASE", "SC_ATAC"])]
        x["uid"] = x.variant_id.str.replace(">", ":", regex=False)
        if keep_uids is not None:
            x = x[x.uid.isin(keep_uids)]
        t = x.groupby(["output_type", "imputation_mode", "track_name"]).size().rename("rows").reset_index()
        t.insert(0, "condition", cond)
        t.insert(0, "dir", d.name)
        tracks.append(t)
        x = x[x.imputation_mode == "observed"]
        g = x.groupby(["uid", "output_type"]).raw_score.mean().reset_index()
        g["readout"] = cond + "|" + g.output_type
        vals.append(g[["uid", "readout", "raw_score"]])
    return vals, tracks


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    v = pd.read_csv(IN / "variants_b2.tsv", sep="\t")
    mv = pd.read_csv(IN / "member_variants_b2.tsv", sep="\t")
    loci = pd.read_csv(IN / "loci_b2.tsv", sep="\t")

    dirs = [(d, d.name.split("_", 1)[1]) for d in sorted(RAW.glob("p1_*"))]
    dirs += [(d, d.name.split("_", 2)[2]) for d in sorted(RAW.glob("p2_b*_*"))]
    vals, tracks = load_dirs(dirs)
    reused = set(v.loc[v.scored_in_stageB, "uid"])
    vb, tb = load_dirs([(d, d.name.split("_", 1)[1]) for d in sorted(RAWB.glob("t?_*"))], keep_uids=reused)
    pd.concat(tracks + tb, ignore_index=True).to_csv(OUT / "stageB2_observed_tracks_used.tsv", sep="\t", index=False)
    s = pd.concat(vals + vb, ignore_index=True).drop_duplicates(["uid", "readout"])
    raw = s.pivot(index="uid", columns="readout", values="raw_score")
    NOTES["variants_with_scores"] = int(raw.shape[0])
    NOTES["reused_from_stageB_with_scores"] = int(len(reused & set(raw.index)))

    lab = [BENCH["BL_ATAC"], BENCH["BL_DNASE"], BENCH["HEP_SCATAC"]]
    have = raw.reindex(columns=lab).notna().all(axis=1)
    vl = v.assign(locus=v.loci.str.split(";")).explode("locus")
    full = vl.groupby("locus").uid.agg(lambda u: bool(have.reindex(u).fillna(False).all()))
    loci["fully_scored"] = loci.locus_id.map(full).fillna(False)
    scored_loci = set(loci.locus_id[loci.fully_scored])
    NOTES["loci_fully_scored_by_tier"] = loci.groupby("tier").fully_scored.sum().astype(int).to_dict()
    NOTES["loci_by_tier"] = loci.groupby("tier").size().astype(int).to_dict()

    bg_all = set(v.loc[v.roles.str.contains("background"), "uid"])
    bgt = pd.read_csv(IN / "background_b2.tsv", sep="\t")
    bg_scored = {u for u, l in zip(bgt.uid, bgt.locus_id) if l in scored_loci} & bg_all
    pct, thr = S.percentiles(raw, bg_scored)
    NOTES["thresholds"] = thr
    sup = (pct >= 0.99).where(pct.notna())

    pv = v.set_index("uid")[["tier", "roles", "loci"]]
    pv = pv[pv.loci.str.split(";").map(lambda ls: any(x in scored_loci for x in ls))]
    pv = pv.join(raw.add_prefix("raw:")).join(pct.add_prefix("pct:"))
    mvi = mv.sort_values("max_pip", ascending=False).drop_duplicates("uid").set_index("uid")
    pv["max_pip"] = pd.to_numeric(mvi.max_pip.reindex(pv.index), errors="coerce")
    pv["member_ancestry"] = mvi.member_ancestry.reindex(pv.index)

    is_mem = pv.roles.str.contains(r"(^|;)member($|;)")
    is_bg = pv.index.isin(bg_scored) & (~is_mem)
    union_df, q99_union, bg_legacy_rate, bg_calib_rate = S.calibrate_union_and_fdr(pv, pct, is_mem, is_bg, "max_pip")
    union_df.to_csv(OUT / "stageB2_union_calibration.tsv", sep="\t", index=False)
    NOTES["union_calibration"] = union_df.to_dict(orient="records")

    pv["sequence_supported"] = pv["calibrated_union_pass"]
    for k, r in BENCH.items():
        pv[f"supported_{k}"] = sup.reindex(pv.index)[r]
    sup_pv = sup.reindex(pv.index)
    pv["supported_readouts"] = [
        ", ".join(f"{k} {pv.at[u, 'pct:' + r]:.3f}" for k, r in BENCH.items() if row[r] == 1.0)
        for u, row in zip(pv.index, sup_pv[list(BENCH.values())].to_dict("records"))
    ]
    pv["supported_conditions"] = sup_pv.apply(
        lambda row: ";".join(
            sorted(
                {
                    c.split("|")[0] + ("(benchmarked)" if c in BENCH.values() else "(not benchmarked)")
                    for c, x in row.items()
                    if x == 1.0
                }
            )
        ),
        axis=1,
    )

    def acc_sign(u):
        z = {np.sign(raw.at[u, r]) for r in lab if sup.at[u, r] == 1.0}
        return "" if not z else ("gain" if z == {1.0} else "loss" if z == {-1.0} else "mixed")

    pv["accessibility_sign_alt_minus_ref"] = [acc_sign(u) if f else "" for u, f in zip(pv.index, pv.sequence_supported)]
    pv["coding"] = "not annotated"
    pv["candidate_not_adopted"] = True
    pv["package_version"] = NOTES["package_version"]
    pv.reset_index().to_csv(OUT / "stageB2_variant_fields.tsv.gz", sep="\t", index=False, compression="gzip")

    bg_share = bg_calib_rate
    NOTES["background_share_calibrated_union"] = {
        "n": int(is_bg.sum()),
        "q99_union": q99_union,
        "bg_calib_rate": bg_calib_rate,
        "bg_legacy_3track_rate": bg_legacy_rate,
    }
    rows = []
    mem = pv[is_mem]
    for key, g in (
        [("all", mem)]
        + [(f"ancestry {a}", x) for a, x in mem.groupby("member_ancestry")]
        + [(f"phase {1 if t == 2 else 2} (tier {t})", x) for t, x in mem.groupby("tier")]
        + [(f"phase {1 if t == 2 else 2} (tier {t}) ancestry {a}", x) for (t, a), x in mem.groupby(["tier", "member_ancestry"])]
    ):
        rows.append(
            {
                "group": key,
                "members": int(len(g)),
                "sequence_supported": int(g.sequence_supported.sum()),
                "prioritized_pip10": int(g.prioritized_pip10.sum()),
                "prioritized_pip50": int(g.prioritized_pip50.sum()),
                "legacy_3track_p99": int(g.legacy_3track_p99.sum()),
                "expected_at_background_share": round(len(g) * bg_share, 1),
                "share": float(g.sequence_supported.mean()) if len(g) else np.nan,
            }
        )
    for k, r in BENCH.items():
        rows.append(
            {
                "group": f"readout {k}",
                "members": int(mem[f"pct:{r}"].notna().sum()),
                "sequence_supported": int((mem[f"pct:{r}"] >= 0.99).sum()),
                "prioritized_pip10": int(((mem[f"pct:{r}"] >= 0.99) & (mem["max_pip"] >= 0.10)).sum()),
                "prioritized_pip50": int(((mem[f"pct:{r}"] >= 0.99) & (mem["max_pip"] >= 0.50)).sum()),
                "legacy_3track_p99": int((mem[f"pct:{r}"] >= 0.99).sum()),
                "expected_at_background_share": round(
                    mem[f"pct:{r}"].notna().sum() * float((pv.loc[is_bg, f"pct:{r}"] >= 0.99).mean()), 1
                ),
            }
        )
    cnt = pd.DataFrame(rows)
    cnt["candidate_not_adopted"] = True
    cnt.to_csv(OUT / "stageB2_counts.tsv", sep="\t", index=False)

    sb_path = L.PKG / "out_20261118/stageB/stageB_variant_fields.tsv.gz"
    sb = pd.read_csv(sb_path, sep="\t", usecols=["uid", "roles", "sequence_supported", "prioritized_pip10", "legacy_3track_p99"])
    sbm = sb[sb.roles.str.contains(r"(^|;)member($|;)")]
    NOTES["combined_with_stageB_separate_line"] = {
        "stageB_members": int(len(sbm)),
        "stageB_sequence_supported": int(sbm.sequence_supported.sum()),
        "stageB_prioritized_pip10": int(sbm.prioritized_pip10.sum()),
        "stageB_legacy_3track": int(sbm.legacy_3track_p99.sum()),
        "stageB2_members": int(len(mem)),
        "stageB2_sequence_supported": int(mem.sequence_supported.sum()),
        "stageB2_prioritized_pip10": int(mem.prioritized_pip10.sum()),
        "stageB2_legacy_3track": int(mem.legacy_3track_p99.sum()),
        "members_in_both_stages": int(len(set(sbm.uid) & set(mem.index))),
    }

    with gzip.open(IN / "credible_set_member_rows_b2.tsv.gz", "rt") as fh:
        cs = pd.read_csv(fh, sep="\t", usecols=["locus_id", "set_id", "method", "pip", "uid", "ancestry"])

    out_loc = []
    for _, lc in loci.sort_values("locus_id").iterrows():
        lid = lc.locus_id
        row = {
            "locus_id": lid,
            "tier": int(lc.tier),
            "genes": lc.genes if isinstance(lc.genes, str) else "",
            "n_credible_sets": int(cs[cs.locus_id == lid].set_id.nunique()),
            "set_ancestries": ",".join(sorted(cs[cs.locus_id == lid].ancestry.unique())),
        }
        if lid not in scored_loci:
            row["sequence_supported_variant"] = "not scored (call cap or locus order)"
            out_loc.append(row)
            continue
        here = mem[mem.loci.str.split(";").map(lambda ls: lid in ls)]
        ss = here[here.sequence_supported]
        row.update(
            {
                "n_members_scored": int(len(here)),
                "sequence_supported_variant": "none (not evidence against the locus)"
                if len(ss) == 0
                else ("one" if len(ss) == 1 else f"several ({len(ss)})"),
                "expected_at_background_share": round(len(here) * bg_share, 2),
                "prioritized_pip10_variants": ";".join(here.index[here.prioritized_pip10]),
                "supported_variants": " | ".join(
                    f"{u} (max PIP {r.max_pip:.3f}; max_bench_pct {r.max_bench_pct:.4f}; "
                    f"{r.supported_conditions}; accessibility {r.accessibility_sign_alt_minus_ref})"
                    for u, r in ss.iterrows()
                ),
            }
        )
        out_loc.append(row)
    lt = pd.DataFrame(out_loc)
    lt["candidate_not_adopted"] = True
    lt[lt.tier == 2].to_csv(OUT / "stageB2_alt_locus_table.tsv", sep="\t", index=False)
    lt[lt.tier == 3].to_csv(OUT / "stageB2_tier3_locus_table.tsv", sep="\t", index=False)

    # B1 (locus-grouped within_set_perm)
    cs2 = cs[cs.uid.isin(pv.index[is_mem]) & cs.locus_id.isin(scored_loci)]
    su = cs2[cs2.method == "susie"].drop_duplicates(["locus_id", "set_id", "uid"]).copy()
    su["set_key"] = su["locus_id"].astype(str) + "||" + su["set_id"].astype(str)
    top = su.sort_values(["set_key", "pip"], ascending=[True, False]).drop_duplicates("set_key")
    su = su.merge(top[["set_key", "uid"]].rename(columns={"uid": "top_uid"}), on="set_key")
    su["block"] = su.top_uid.map(S.block)
    pct_m = pct.reindex(columns=list(dict.fromkeys(lab + S.CTRL_SCATAC + [c for c in pct.columns])))
    b1 = S.within_set_perm_locus_grouped(su[["locus_id", "set_id", "uid", "pip", "block"]], pct_m, L.SEED)
    b1["candidate_not_adopted"] = True
    b1.to_csv(OUT / "stageB2_B1.tsv", sep="\t", index=False)

    # B2
    seqsup_calib = pv["calibrated_union_pass"].to_dict()
    seqsup_pip10 = pv["prioritized_pip10"].to_dict()
    seqsup_legacy = pv["legacy_3track_p99"].to_dict()
    lm = {lid: list(set(g.uid)) for lid, g in cs2.groupby("locus_id")}
    lb = {lid: list(set(g.uid)) for lid, g in bgt[bgt.locus_id.isin(scored_loci)].groupby("locus_id")}
    r2 = S.b2(lm, lb, seqsup_calib)
    r2["calibrated_union_pip10"] = S.b2(lm, lb, seqsup_pip10)
    r2["legacy_uncalibrated_3track"] = S.b2(lm, lb, seqsup_legacy)
    r2["candidate_not_adopted"] = True
    (OUT / "stageB2_B2.json").write_text(json.dumps(r2, indent=1))

    # B3 input (np.nan when Currin entry missing)
    C = S.CURRIN
    lead = pd.read_csv(
        C / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
        sep="\t",
        usecols=["lead_variant_ID", "peak_ID", "beta"],
    )
    prox = pd.read_csv(
        C / "caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz",
        sep="\t",
        usecols=["proxy_ID", "lead_ID", "caPeak", "proxy_caQTL_beta"],
    )
    chk = prox[prox.proxy_ID == prox.lead_ID].merge(
        lead, left_on=["lead_ID", "caPeak"], right_on=["lead_variant_ID", "peak_ID"]
    )
    agree = float((np.sign(chk.proxy_caQTL_beta) == np.sign(chk.beta)).mean())
    eff = lead.rename(columns={"lead_variant_ID": "id", "beta": "b"})[["id", "b"]]
    if agree >= 0.99:
        eff = pd.concat([eff, prox.rename(columns={"proxy_ID": "id", "proxy_caQTL_beta": "b"})[["id", "b"]]])
    eff = eff.assign(a=eff.b.abs()).sort_values("a", ascending=False).drop_duplicates("id")
    ca = dict(zip(eff.id, eff.b))
    NOTES["B3a_proxy_check"] = {"lead_rows": int(len(chk)), "agreement": agree, "proxies_used": bool(agree >= 0.99)}
    d1 = pd.read_csv(S.D1, sep="\t")
    d1 = d1[(d1.scope == "POOLED_30_DONORS") & (d1.lineage == "Hepatocyte")]
    ai = {u: float(abs(m / e) >= 1.96) for u, m, e in zip(d1.uid, d1.mean_log2_alt_over_ref, d1.se_log2) if e > 0}
    ai.update({flip(u): x for u, x in list(ai.items()) if flip(u) not in ai})
    mp = pd.read_csv(S.MPRA, sep="\t", dtype=str)
    mp = mp[mp.context.isin(["HepG2_ctrl", "HepG2_PAOA"]) & (mp.mpra_barcode_mapped.str.upper() == "TRUE")]
    mp["uid"] = mp.chr + ":" + mp.pos_hg38 + ":" + mp.ref_hg38 + ":" + mp.alt_hg38
    mpd = mp.assign(d=mp.dav.astype(float)).groupby("uid").d.max().to_dict()
    mpd.update({flip(u): x for u, x in list(mpd.items()) if flip(u) not in mpd})

    def ca_out(u):
        b = ca.get(u, -ca[flip(u)] if flip(u) in ca else None)
        if b is None:
            return np.nan
        val = raw.at[u, BENCH["BL_ATAC"]] if u in raw.index else np.nan
        return np.nan if not np.isfinite(val) else float(np.sign(b) == np.sign(val))

    arms = {
        "B3a_caqtl_matching_sign": (sup[BENCH["BL_ATAC"]], ca_out),
        "B3b_hepatocyte_allelic_imbalance": (sup[BENCH["HEP_SCATAC"]], lambda u: ai.get(u, np.nan)),
        "B3c_mpra_allelic_activity": (pd.Series({u: float(x) for u, x in seqsup_calib.items()}), lambda u: mpd.get(u, np.nan)),
    }
    mem_cs = cs2[["locus_id", "set_id", "uid"]].drop_duplicates()
    b3 = []
    for arm, (sp, outcome) in arms.items():
        d = mem_cs.assign(supported=mem_cs.uid.map(sp)).dropna(subset=["supported"])
        ok = d.groupby("set_id").supported.agg(lambda x: x.min() == 0 and x.max() == 1)
        d = d[d.set_id.isin(ok[ok].index)].drop_duplicates(["locus_id", "uid"])
        d["outcome"] = d.uid.map(outcome)
        d = d.dropna(subset=["outcome"]).assign(arm=arm)
        b3.append(d[["arm", "locus_id", "uid", "supported", "outcome"]])
    pd.concat(b3).to_csv(OUT / "stageB2_B3_input.tsv", sep="\t", index=False)
    NOTES["B3a_members_with_caqtl_entry"] = int(sum((u in ca) or (flip(u) in ca) for u in mem_cs.uid.unique()))
    NOTES["B3b_members_tested"] = int(sum(u in ai for u in mem_cs.uid.unique()))
    NOTES["B3c_members_tested"] = int(sum(u in mpd for u in mem_cs.uid.unique()))

    ns = loci[~loci.fully_scored]
    NOTES["loci_not_scored"] = {
        "n": int(len(ns)),
        "by_tier": ns.groupby("tier").size().astype(int).to_dict(),
        "member_variants": int(mv[mv.primary_locus.isin(ns.locus_id)].shape[0]),
    }
    (OUT / "stageB2_notes.json").write_text(json.dumps(NOTES, indent=1, default=str))
    L.log(f"Stage B2 20261118 complete in {OUT}")


if __name__ == "__main__":
    main()
