#!/usr/bin/env python3
"""Model A-final Gate A design export for the independent design simulation (Codex C5/C8).

Outcome-free: uses total depth n = REF + ALT only (no allele fraction is formed), stage and
covariates. Written outside the restricted directory, so two things are coarsened:
  - the genotype prior of each row uses the cohort's 1kGP superpopulation mixture (headcount of
    the set's ancestry calls), not the individual's own ancestry call;
  - ancestry axes are not exported per individual; their effect on the within-cohort R2 of S is
    given per cohort (with and without the axes) in design_cohort_moments.json.

Outputs in --out:
  design_rows.tsv.gz        one row per set x individual x gene: the chosen tag (spec v1 1.2: largest
                            n among the gene's kept tags, ties to the smallest position; FFPE transition
                            drop applied), its s_t, n, and the merged prior weights P_H, P_R, P_A (2.1)
  design_individuals.tsv    every unit individual of the set (with or without rows): cohort, cluster_id
                            (first-degree / same-individual components), S, z columns (spec 3.3 without
                            the axes: b = log max(e_model, 1e-4), d = dup1_sub_rpk2to20, ffpe =
                            ctga_excess_nontarget_clean where reliable else cohort mean, sex, age/10; v1.4), and the stage-law
                            depth predictors log1p(candidate rows), log1p(median candidate n)
  design_cohort_moments.json  per set x cohort: individuals, R2 of S on z (without and with the two axes),
                            and the rule R_c
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

STAGED5 = ["GSE213621", "GSE135251", "GSE162694", "GSE130970", "GSE240729"]
JAPANESE = {"GSE167523", "GSE174478", "GSE193066"}
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS", "JPT"]


def components(pairs, members):
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in pairs:
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    return {m: find(m) for m in members}


def r2(y, X):
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    y, X = y[ok], X[ok]
    if len(y) < X.shape[1] + 2 or np.var(y) == 0:
        return float("nan")
    X1 = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return float(1 - np.var(y - X1 @ beta) / np.var(y))


def main():
    ap = argparse.ArgumentParser()
    for k in ["crosswalk", "ad", "tag-table", "tag-freq", "ancestry", "library-qc", "dup", "metadata",
              "kinship", "rules", "out"]:
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    xw = pd.read_csv(a.crosswalk, sep="\t", dtype=str)
    xw = xw[xw["unit_library"] == "True"]
    staged = (xw["C"] == "0") & xw["S"].notna() & (xw["S"] != "")
    dev = (xw["run_role"] == "development") & xw["cohort"].isin(STAGED5) & staged
    sea = (xw["run_role"] == "sealed") & staged & xw["cohort"].isin(STAGED5 + ["GSE174478"])
    units = pd.concat([xw[dev].assign(set="development"), xw[sea].assign(set="sealed")])[
        ["run", "cohort", "individual_id", "set", "S"]]

    tags = pd.read_csv(a.tag_table, sep="\t", dtype={"pos": int})
    drop_cols = [c for c in tags.columns if c.startswith("drop_in_")]
    tags = tags[tags["kept_primary"]]
    ad = pd.read_csv(a.ad, sep="\t", usecols=["run", "chrom", "pos", "ref_reads", "alt_reads"])
    ad["n"] = ad["ref_reads"] + ad["alt_reads"]
    ad = ad.drop(columns=["ref_reads", "alt_reads"]).merge(units, on="run")
    rows = ad.merge(tags[["gene_id", "chrom", "pos", "s_t", "transition"] + drop_cols], on=["chrom", "pos"])
    for c in drop_cols:
        coh = c[len("drop_in_"):]
        rows = rows[~((rows["cohort"] == coh) & rows[c])]
    rows = rows.sort_values(["run", "gene_id", "n", "pos"], ascending=[True, True, False, True])
    chosen = rows.drop_duplicates(["run", "gene_id"]).copy()

    # cohort-mixture prior (set x cohort headcount of ancestry calls)
    anc = pd.read_csv(a.ancestry, sep="\t", usecols=["run", "ancestry_call", "nearest_superpop", "intermediate", "axis_EUR_AFR", "axis_EUR_EAS"])
    anc["pop"] = np.where(anc["intermediate"].astype(str) == "True", anc["nearest_superpop"], anc["ancestry_call"])
    u2 = units.merge(anc, on="run", how="left")
    u2.loc[u2["cohort"].isin(JAPANESE) & (u2["pop"] == "EAS"), "pop"] = "JPT"
    mix = u2.dropna(subset=["pop"]).groupby(["set", "cohort"])["pop"].value_counts(normalize=True).rename("w").reset_index()
    fr = pd.read_csv(a.tag_freq, sep="\t", dtype={"pos": int})
    parts = []
    for (s_, c), m in mix.groupby(["set", "cohort"]):
        ph = pr = pa = 0.0
        for p, w in zip(m["pop"], m["w"]):
            pref, palt = fr[f"p_ref_{p}"], fr[f"p_alt_{p}"]
            px = fr[f"other_alt_freq_{p}"] + fr[f"deleted_freq_{p}"] + fr[f"conflict_freq_{p}"]
            n_hap = fr[f"n_hap_{p}"].clip(lower=2)
            lo = 1 / n_hap
            pref, palt = pref.clip(lo, 1 - lo), palt.clip(lo, 1 - lo)
            ph = ph + w * 2 * pref * palt
            pr = pr + w * (pref ** 2 + 2 * pref * px)
            pa = pa + w * (palt ** 2 + 2 * palt * px)
        parts.append(pd.DataFrame({"set": s_, "cohort": c, "chrom": fr["chrom"], "pos": fr["pos"], "P_H": ph, "P_R": pr, "P_A": pa}))
    prior = pd.concat(parts, ignore_index=True)
    chosen = chosen.merge(prior, on=["set", "cohort", "chrom", "pos"], how="left")
    cols = ["set", "cohort", "individual_id", "gene_id", "chrom", "pos", "s_t", "n", "P_H", "P_R", "P_A"]
    chosen[cols].sort_values(cols[:4]).to_csv(out / "design_rows.tsv.gz", sep="\t", index=False, compression="gzip")

    # individuals and z (spec 3.3 without the axes)
    # spec v1.4 z binding: b = log max(e_model, 1e-4) (B-QC r1, e_model per 2.3); FFPE per library only where
    # its split-half reliability allows it (z_ffpe_per_library_allowed), else the cohort mean; c_i dropped
    qc = pd.read_csv(a.library_qc, sep="\t", usecols=["run", "e_model", "expr_sex", "ctga_excess_nontarget_clean",
                                                     "z_ffpe_per_library_allowed"])
    dup = pd.read_csv(a.dup, sep="\t", usecols=["run", "dup1_sub_rpk2to20"])
    md = pd.read_csv(a.metadata, usecols=["sample_id", "sex", "age"]).rename(columns={"sample_id": "run"})
    ind = u2.merge(qc, on="run", how="left").merge(dup, on="run", how="left").merge(md, on="run", how="left")
    sex = ind["sex"].astype(str).str.upper().str[0].where(ind["sex"].notna())
    sex = sex.where(sex.isin(["M", "F"]), ind["expr_sex"].where(ind["expr_sex"].isin(["M", "F"])))
    ind["z_sex_female"] = (sex == "F").astype(float).where(sex.notna())
    ind["z_b_log_e"] = np.log(ind["e_model"].clip(lower=1e-4))
    ind["z_d_dup"] = ind["dup1_sub_rpk2to20"]
    allowed = ind["z_ffpe_per_library_allowed"].astype(str) == "True"
    cohort_mean = ind.groupby(["set", "cohort"])["ctga_excess_nontarget_clean"].transform("mean")
    ind["z_ffpe"] = np.where(allowed, ind["ctga_excess_nontarget_clean"], cohort_mean)
    ind["z_age10"] = pd.to_numeric(ind["age"], errors="coerce") / 10
    cand = chosen.groupby("run")["n"].agg(["size", "median"])
    ind = ind.merge(cand, left_on="run", right_index=True, how="left")
    ind["stage_law_log1p_candidate_rows"] = np.log1p(ind["size"].fillna(0))
    ind["stage_law_log1p_median_n"] = np.log1p(ind["median"].fillna(0))
    kin = pd.read_csv(a.kinship, sep="\t")
    kin = kin[kin["relation"].isin(["first_degree", "same_individual"])]
    run2ind = dict(zip(xw["run"], xw["individual_id"]))
    pairs = [(run2ind.get(x), run2ind.get(y)) for x, y in zip(kin["#IID1"], kin["IID2"])]
    comp = components([p for p in pairs if None not in p], ind["individual_id"].tolist())
    ind["cluster_id"] = ind["individual_id"].map(comp)
    zc = ["z_b_log_e", "z_d_dup", "z_ffpe", "z_sex_female", "z_age10"]
    keep = ["set", "cohort", "individual_id", "cluster_id", "S"] + zc + ["stage_law_log1p_candidate_rows", "stage_law_log1p_median_n"]
    ind[keep].sort_values(["set", "cohort", "individual_id"]).to_csv(out / "design_individuals.tsv", sep="\t", index=False)

    rules = pd.read_csv(a.rules, sep="\t")
    mom = {}
    for (s_, c), g in ind.groupby(["set", "cohort"]):
        S = pd.to_numeric(g["S"]).to_numpy(float)
        Z = g[zc].apply(lambda v: v.fillna(v.mean())).to_numpy(float)
        Z = Z[:, np.nanstd(Z, axis=0) > 0] if len(Z) else Z
        A = g[["axis_EUR_AFR", "axis_EUR_EAS"]].apply(lambda v: v.fillna(v.mean())).to_numpy(float)
        r = rules[(rules["cohort"] == c) & (rules["set"] == s_)].iloc[0]
        mom[f"{s_}|{c}"] = {"individuals": int(len(g)), "individuals_with_rows": int(g["size"].notna().sum()),
                            "r2_S_on_z_without_axes": r2(S, Z), "r2_S_on_z_with_axes": r2(S, np.column_stack([Z, A])),
                            "rule": [int(r["min_dp"]), int(r["min_minor_reads"]), int(r["min_minor_pct"])],
                            "clusters_with_more_than_one": int((g.groupby("cluster_id").size() > 1).sum()),
                            "missing_share": {k: float(g[k].isna().mean()) for k in zc}}
    xc = ind.groupby("cluster_id")["cohort"].nunique()
    mom["clusters_crossing_cohorts"] = int((xc > 1).sum())
    (out / "design_cohort_moments.json").write_text(json.dumps(mom, indent=1))
    print(json.dumps(mom, indent=1))


if __name__ == "__main__":
    main()
