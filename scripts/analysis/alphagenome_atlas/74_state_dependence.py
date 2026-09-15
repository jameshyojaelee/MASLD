#!/usr/bin/env python3
"""Step 74 (P3): does disease state modify allelic accessibility at Resource variants?

An allelic effect that differs between MASH and normal donors is a state-conditional regulatory variant, and
a static sequence model cannot represent one. That is the MASLD-specific question this substrate can ask of
the Atlas's predictions.

Per-site testing is abandoned as untestable at this n rather than reported as a null: with >= 3 het donors
per group a per-site rank test has a minimum two-sided p of 0.1. The powered question is global -- whether
any state-dependence exists beyond noise -- with sites as repeated measures and donors as the unit. A global
signal would not license naming any variant.

Prespecification: 74_state_dependence_prespec.json (filters, seeds, null, written predictions).
usage: 74_state_dependence.py [gse281367|gse244832]
Outputs (tables/): <prefix>state_dependence.json, <prefix>state_deltas.tsv
"""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import sys

import numpy as np
import pandas as pd

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
PRESPEC = la.SCRIPT_DIR / "74_state_dependence_prespec.json"
CONTROL_LABELS = {"NORMAL", "HEALTHY", "CONTROL"}
COHORTS = {
    "gse281367": {"counts": "allelic_donor_counts.tsv.gz", "prefix": "",
                  "meta": "Analysis/ATAC/Human_External/pseudobulk/hep_pseudobulk_coldata_GSE281367.tsv"},
    "gse244832": {"counts": "gse244832_allelic_donor_counts.tsv.gz", "prefix": "gse244832_",
                  "meta": "Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv"},
}
DRAWS = 2000
SEED = 20260909
MIN_PER_GROUP = 2
DEPTH_GAP_LIMIT = 0.15


def site_deltas(rows: pd.DataFrame, min_per_group: int = MIN_PER_GROUP) -> pd.DataFrame:
    """Per site: mean donor log2 ratio in disease minus mean in control, donors as the unit."""
    g = rows.groupby(["uid", "grp"])["log2"].agg(["mean", "size"]).unstack()
    if ("size", "disease") not in g or ("size", "control") not in g:
        return pd.DataFrame(columns=["uid", "delta", "n_disease", "n_control"])
    ok = (g[("size", "disease")].fillna(0) >= min_per_group) & (g[("size", "control")].fillna(0) >= min_per_group)
    g = g[ok]
    return pd.DataFrame({"uid": g.index,
                         "delta": (g[("mean", "disease")] - g[("mean", "control")]).values,
                         "n_disease": g[("size", "disease")].values,
                         "n_control": g[("size", "control")].values})


def permuted_group_maps(donors: list, base: dict, draws: int, seed: int, strata: dict | None = None) -> list:
    """Donor-to-group assignments with group sizes preserved, optionally shuffled within strata.

    Shuffling donor labels keeps every within-donor and between-site correlation intact; shuffling sites
    would destroy the very structure the statistic lives in.
    """
    rng = np.random.default_rng(seed)
    donors = list(donors)
    out = []
    if strata is None:
        labels = [base[d] for d in donors]
        for _ in range(draws):
            out.append(dict(zip(donors, rng.permutation(labels))))
        return out
    by_stratum = {}
    for d in donors:
        by_stratum.setdefault(strata[d], []).append(d)
    for _ in range(draws):
        mp = {}
        for members in by_stratum.values():
            labels = [base[d] for d in members]
            mp.update(dict(zip(members, rng.permutation(labels))))
        out.append(mp)
    return out


def group_depth_gap(rows: pd.DataFrame) -> float:
    """|mean log10 depth (disease) - mean log10 depth (control)| over donors.

    The +0.5 pseudocount shrinks a log ratio toward zero at low depth, so a depth-group association would
    manufacture a delta on its own.
    """
    per = rows.groupby(["donor", "grp"])["depth"].mean().reset_index()
    m = per.groupby("grp")["depth"].apply(lambda s: float(np.mean(np.log10(s))))
    if "disease" not in m or "control" not in m:
        return float("nan")
    return float(abs(m["disease"] - m["control"]))


def main() -> None:
    name = (sys.argv[1] if len(sys.argv) > 1 else "gse281367").lower()
    if name not in COHORTS:
        raise la.ContractError(f"unknown cohort {name!r}")
    cfg = COHORTS[name]
    spec = json.loads(PRESPEC.read_text())
    digest = hashlib.sha256(PRESPEC.read_bytes()).hexdigest()
    shutil.copy(PRESPEC, ROOT / "74_state_dependence_prespec.json")

    counts = pd.read_csv(TABLES / cfg["counts"], sep="\t")
    meta = pd.read_csv(la.PROJECT / cfg["meta"], sep="\t").rename(columns={"donor_id": "donor"})
    meta = meta[["donor", "condition"]].drop_duplicates()
    d = counts[counts["is_het"]].merge(meta, on="donor", how="left")
    if d["condition"].isna().any():
        raise la.ContractError(f"donors without a condition label: {sorted(d.loc[d.condition.isna(), 'donor'].unique())}")
    d["grp"] = np.where(d["condition"].astype(str).str.upper().isin(CONTROL_LABELS), "control", "disease")
    d["log2"] = np.log2((d["n_alt"] + 0.5) / (d["n_ref"] + 0.5))
    d["depth"] = d["n_ref"] + d["n_alt"]

    obs = site_deltas(d)
    summary = {"cohort": name, "prespec_sha256": digest, "n_donors": int(d["donor"].nunique()),
               "group_sizes": d.groupby("grp")["donor"].nunique().to_dict(),
               "n_eligible_sites": int(len(obs)),
               "per_site_testing": spec["power_statement_written_before_results"]["per_site_testing"]}

    if len(obs) < 50:
        summary["verdict"] = "cannot conclude: fewer than 50 eligible sites"
        json.dump(summary, (TABLES / f"{cfg['prefix']}state_dependence.json").open("w"), indent=1, default=float)
        la.log(f"P3 state [{name}]: {len(obs)} eligible sites, cannot conclude")
        return

    gap = group_depth_gap(d)
    summary["group_depth_gap_log10"] = gap
    strata = None
    if gap == gap and gap > DEPTH_GAP_LIMIT:
        per = d.groupby("donor")["depth"].mean()
        strata = pd.qcut(per, q=3, labels=False, duplicates="drop").to_dict()
        summary["permutation_stratified_on_depth_tertile"] = True
    else:
        summary["permutation_stratified_on_depth_tertile"] = False

    base = d.drop_duplicates("donor").set_index("donor")["grp"].to_dict()
    donors = sorted(base)
    obs_mean_abs = float(np.mean(np.abs(obs["delta"])))
    obs_sd = float(np.std(obs["delta"], ddof=1))
    null_mean_abs, null_sd = [], []
    for mp in permuted_group_maps(donors, base, DRAWS, SEED, strata):
        perm = d.copy()
        perm["grp"] = perm["donor"].map(mp)
        pd_ = site_deltas(perm)
        if len(pd_) == 0:
            continue
        null_mean_abs.append(float(np.mean(np.abs(pd_["delta"]))))
        null_sd.append(float(np.std(pd_["delta"], ddof=1)))
    null_mean_abs = np.asarray(null_mean_abs)
    null_sd = np.asarray(null_sd)

    summary.update({
        "observed_mean_abs_delta": obs_mean_abs,
        "null_mean_abs_delta_mean": float(null_mean_abs.mean()),
        "null_mean_abs_delta_ci95": [float(np.percentile(null_mean_abs, 2.5)), float(np.percentile(null_mean_abs, 97.5))],
        "exceedance_count_mean_abs": int(np.sum(null_mean_abs >= obs_mean_abs)),
        "p_mean_abs": float((np.sum(null_mean_abs >= obs_mean_abs) + 1) / (len(null_mean_abs) + 1)),
        "observed_sd_delta": obs_sd,
        "null_sd_delta_mean": float(null_sd.mean()),
        "exceedance_count_sd": int(np.sum(null_sd >= obs_sd)),
        "p_sd": float((np.sum(null_sd >= obs_sd) + 1) / (len(null_sd) + 1)),
        "n_draws": int(len(null_mean_abs)),
    })
    summary["verdict"] = ("global state-dependence beyond the donor-label null"
                          if summary["p_mean_abs"] < 0.05 or summary["p_sd"] < 0.05
                          else "no global state-dependence detectable at this n")
    summary["claim_boundary"] = spec["claim_boundary"]
    obs.to_csv(TABLES / f"{cfg['prefix']}state_deltas.tsv", sep="\t", index=False)
    json.dump(summary, (TABLES / f"{cfg['prefix']}state_dependence.json").open("w"), indent=1, default=float)
    la.log(f"P3 state [{name}]: {len(obs)} sites, mean|delta| {obs_mean_abs:.4f} vs null {null_mean_abs.mean():.4f}, "
           f"exceedance {summary['exceedance_count_mean_abs']}/{len(null_mean_abs)}")


if __name__ == "__main__":
    main()
