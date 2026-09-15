#!/usr/bin/env python3
"""Step 44 (P4): the prespecified tests on the saturation features.

Prespecification 43_saturation_prespec.json, written before any feature existed. Families:
  F1  are the H3K27ac regions the released chromatin predictor can predict more sequence-sensitive than the
      ones it cannot? Matched draws within GC x width x promoter x signal-mean x signal-sd strata.
  F2  is sensitivity more localised in promoter regions than elsewhere? Same matching machinery.
  F3  which TF tracks carry the sensitivity, and does the composition differ between the two groups?

The primary scale is the relative raw effect; every F1/F2 statistic is also reported on the saturating
quantile scale so a reader can see the difference the scale makes.

Outputs (tables/): saturation_tests.json, saturation_tf_composition.tsv
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la

PRESPEC = pathlib.Path(__file__).resolve().parent / "43_saturation_prespec.json"
COVARIATES = la.track0_root() / "tables" / "region_covariates.tsv.gz"
N_DRAWS = 10_000
SEED = 20260914
GROUPS_PRIMARY = ("liver_h3k27ac", "liver_atac", "liver_dnase")


def bins(values, n: int) -> np.ndarray:
    """Quantile bins that tolerate ties and missing values (missing gets its own bin, -1)."""
    v = pd.to_numeric(pd.Series(values), errors="coerce")
    out = pd.Series(-1, index=v.index, dtype=int)
    ok = v.notna()
    if ok.sum():
        ranks = v[ok].rank(method="average", pct=True)
        out.loc[ok] = np.clip((ranks * n).astype(int), 0, n - 1)
    return out.to_numpy()


def stratum_labels(df: pd.DataFrame, use_promoter: bool = True) -> np.ndarray:
    """One label per region: GC decile, width decile, signal mean and sd quintiles, and promoter status.

    `use_promoter=False` is for the contrast whose GROUPS are promoter and non-promoter regions: matching on
    a variable that defines the groups would leave no control for anything.
    """
    parts = [bins(df["gc"], 10), bins(df["width"], 10),
             bins(df["signal_mean"], 5), bins(df["signal_sd"], 5)]
    if use_promoter:
        parts.append(pd.to_numeric(df["promoter"], errors="coerce").fillna(-1).astype(int).to_numpy())
    return np.array(["|".join(str(p[i]) for p in parts) for i in range(len(df))])


def matched_draw_null(values, strata, in_group, n_draws: int, seed: int, statistic: str = "median") -> dict:
    """Median difference (group - matched control) against draws that keep the group's stratum composition.

    Each draw replaces every group region with a random region of the SAME stratum from outside the group,
    so the null holds GC, width, promoter status and measured signal fixed. Regions in a stratum with no
    outside region are dropped from both the statistic and the null, and counted.
    """
    v = np.asarray(values, dtype=float)
    s = np.asarray(strata)
    g = np.asarray(in_group, dtype=bool)
    pool = {}
    for lab in np.unique(s[~g]):
        idx = np.where((~g) & (s == lab))[0]
        if idx.size:
            pool[lab] = idx
    usable = g & np.isfinite(v) & np.array([lab in pool for lab in s])
    if usable.sum() == 0:
        return {"n_group": 0, "n_dropped_no_match": int(g.sum()), "observed": None}
    agg = {"median": np.nanmedian, "mean": np.nanmean}[statistic]
    obs_group = float(agg(v[usable]))
    rng = np.random.default_rng(seed)
    labs = s[usable]
    flat = np.concatenate([pool[lab] for lab in sorted(pool)])
    starts, sizes = {}, {}
    at = 0
    for lab in sorted(pool):
        starts[lab], sizes[lab] = at, pool[lab].size
        at += pool[lab].size
    st = np.array([starts[lab] for lab in labs])
    sz = np.array([sizes[lab] for lab in labs])
    draws = np.empty(n_draws)
    for d in range(n_draws):
        picks = flat[st + (rng.random(st.size) * sz).astype(int)]
        draws[d] = float(agg(v[picks]))
    diff = obs_group - draws
    return {"n_group": int(usable.sum()), "n_dropped_no_match": int(g.sum() - usable.sum()),
            "statistic": statistic, "observed_group": obs_group, "null_median": float(np.median(draws)),
            "null_lo": float(np.quantile(draws, 0.025)), "null_hi": float(np.quantile(draws, 0.975)),
            "observed_minus_null": float(np.median(diff)),
            "n_draws_at_least_as_extreme": int((np.abs(draws - np.median(draws)) >= abs(obs_group - np.median(draws))).sum()),
            "n_draws": n_draws,
            "exceedance_p": float(((np.abs(draws - np.median(draws)) >= abs(obs_group - np.median(draws))).sum() + 1) / (n_draws + 1))}


def cmh_test(present, in_group, strata) -> dict:
    """Cochran-Mantel-Haenszel test of a binary indicator between two groups across strata."""
    from scipy.stats import chi2
    x = np.asarray(present, dtype=float)
    g = np.asarray(in_group, dtype=bool)
    s = np.asarray(strata)
    num, var = 0.0, 0.0
    for lab in np.unique(s):
        m = s == lab
        n = m.sum()
        n1 = (m & g).sum()
        if n1 == 0 or n1 == n or n < 2:
            continue
        t = x[m].sum()
        a = x[m & g].sum()
        num += a - n1 * t / n
        var += n1 * (n - n1) * t * (n - t) / (n * n * (n - 1))
    if var <= 0:
        return {"chi2": None, "p": 1.0}
    stat = max(abs(num) - 0.5, 0.0) ** 2 / var
    return {"chi2": float(stat), "p": float(chi2.sf(stat, 1)), "excess_in_group": float(num)}

def spearman_within_strata(x, y, strata) -> dict:
    """Spearman correlation of x and y after ranking within each stratum (strata of < 3 regions skipped)."""
    df = pd.DataFrame({"x": pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(),
                       "y": pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(), "s": np.asarray(strata)})
    df = df.dropna()
    df = df[df.groupby("s")["x"].transform("size") >= 3]
    if len(df) < 3:
        return {"n": int(len(df)), "rho": None}
    rx = df.groupby("s")["x"].rank(pct=True)
    ry = df.groupby("s")["y"].rank(pct=True)
    return {"n": int(len(df)), "n_strata": int(df["s"].nunique()),
            "rho": float(np.corrcoef(rx, ry)[0, 1])}


def latest_features() -> pathlib.Path:
    """The newest complete feature run; pilots and superseded runs never count."""
    root = la.out_root().parent
    runs = [p for p in sorted(root.glob("p4-features-*/tables/saturation_region_features.tsv.gz"))
            if "PILOT" not in p.parts[-3].upper() and not (p.parents[1] / "SUPERSEDED.txt").exists()]
    if not runs:
        raise la.ContractError("no complete p4-features run")
    return runs[-1]


def main() -> None:
    import hashlib
    tables = la.out_root() / "tables"
    fpath = latest_features()
    feats = pd.read_csv(fpath, sep="\t", low_memory=False)
    cov = pd.read_csv(COVARIATES, sep="\t", low_memory=False)
    u1 = feats[feats["universe"] == "U1_h3k27ac"]
    df = u1.merge(cov, on="region_key", how="left", suffixes=("", "_cov"), validate="one_to_one")
    out = {"prespec_sha256": hashlib.sha256(PRESPEC.read_bytes()).hexdigest(), "features": str(fpath),
           "n_u1_regions": int(len(u1)), "archive_states": df["archive_state"].value_counts().to_dict(),
           "n_without_covariates": int(df["gc"].isna().sum())}
    df = df[(df["archive_state"] == "ok") & df["gc"].notna()].reset_index(drop=True)
    group = df["reliable_shipped_form"].astype(str).eq("True").to_numpy()
    strata = stratum_labels(df)
    out["n_analysed"], out["n_predictable"] = int(len(df)), int(group.sum())
    f1 = {}
    for grp in GROUPS_PRIMARY:
        for scale in ("rel", "q"):
            col = f"{grp}_{scale}_mean"
            f1[col] = matched_draw_null(df[col], strata, group, N_DRAWS, SEED)
            f1[col]["skill_vs_sensitivity_within_strata"] = spearman_within_strata(
                df["skill_shipped_mean"], df[col], strata)
    out["F1_predictability"] = f1
    prom = pd.to_numeric(df["promoter"], errors="coerce").eq(1).to_numpy()
    s_np = stratum_labels(df, use_promoter=False)
    out["F2_positional"] = {f"{g}_{sc}_share_best50": matched_draw_null(df[f"{g}_{sc}_share_best50"], s_np, prom, N_DRAWS, SEED)
                            for g in GROUPS_PRIMARY for sc in ("rel", "q")}
    out["P5_signal_mean"] = {f"{g}_rel_mean": spearman_within_strata(df["signal_mean"], df[f"{g}_rel_mean"], np.zeros(len(df)))
                             for g in GROUPS_PRIMARY}
    tf = pd.read_csv(fpath.parent / "saturation_region_tf_top5.tsv.gz", sep="\t", low_memory=False)
    keep = set(df["region_key"])
    tf = tf[tf["region_key"].isin(keep)]
    idx = {k: i for i, k in enumerate(df["region_key"])}
    rows, ps = [], []
    for (grp, track), sub in tf.groupby(["group", "track_name"]):
        present = np.zeros(len(df))
        present[[idx[k] for k in sub["region_key"]]] = 1.0
        if present.sum() < 50:
            continue
        res = matched_draw_null(present, strata, group, 1000, SEED, statistic="mean")
        cmh = cmh_test(present, group, strata)
        rows.append({"group": grp, "track_name": track, "transcription_factor": sub["transcription_factor"].iloc[0],
                     "share_all_regions": float(present.mean()), "share_predictable": res["observed_group"],
                     "share_matched_null": res["null_median"], "cmh_p": cmh["p"]})
        ps.append(cmh["p"])
    if rows:
        from scipy import stats as sps
        qs = sps.false_discovery_control(ps, method="bh")
        for r, q in zip(rows, qs):
            r["bh_q"] = float(q)
    la.write_tsv_once(tables / "saturation_tf_composition.tsv", rows,
                      ["group", "track_name", "transcription_factor", "share_all_regions", "share_predictable",
                       "share_matched_null", "cmh_p", "bh_q"])
    out["F3_tf_composition"] = {"n_tracks_tested": len(rows),
                                "n_bh_significant": sum(1 for r in rows if r.get("bh_q", 1) < 0.05),
                                "top_tracks_overall": [r["track_name"] for r in
                                                       sorted(rows, key=lambda r: -r["share_all_regions"])[:10]]}
    json.dump(out, (tables / "saturation_tests.json").open("w"), indent=1, default=float)
    la.log(f"step 44: {len(df)} regions, {int(group.sum())} predictable; F3 {len(rows)} TF tracks")


if __name__ == "__main__":
    main()
