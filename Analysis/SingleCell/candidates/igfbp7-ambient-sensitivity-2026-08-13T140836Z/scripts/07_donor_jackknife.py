#!/usr/bin/env python
"""Leave-one-donor-out check on the ambient SENSITIVITY result.

WHY: the frozen primary donor model has n=64 and max leverage 0.547, far above
the 6/64 = 0.094 average. If a single donor carried the module-8 beta collapse
the verdict would be an outlier story, not an ambient story. This refits the
frozen model dropping one donor at a time and reports the range of
%-stage-beta-retained for the focus programs.

SCOPE (binding): SENSITIVITY ANALYSIS ONLY. Reads this candidate's own step-03
donor scores plus the read-only donor metadata. Refits nothing frozen and
rediscovers no program.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(os.environ["MASLD_PROJECT_ROOT"])
CAND = Path(os.environ["CAND_ROOT"])
RES = CAND / "results"

META = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
PAIRINGS = {
    "GSE244832": ROOT / "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": ROOT / "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": ROOT / "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": ROOT / "data/GSE136103/metadata/donor_pairing.csv",
}
PRIMARY_STAGES = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}
FOCUS = {8: "Stromal ECM (IGFBP7)", 20: "Ductular injury (BICC1)",
         2: "Sinusoidal endothelial (STAB2)", 1: "Leukocyte / immune (DOCK2)",
         12: "Secretory plasma protein (ALB)", 21: "Xenobiotic / drug metab (CYP)",
         16: "PPARa lipid metabolism"}


def donor_map() -> dict[str, str]:
    out: dict[str, str] = {}
    for dataset, path in PAIRINGS.items():
        table = pd.read_csv(path, dtype=str)
        for row in table.itertuples(index=False):
            donor = f"{dataset}_{row.donor_id}"
            for sample in str(row.rna_srrs).split(";"):
                if sample.strip():
                    out[sample.strip()] = donor
    return out


def donor_metadata(mapping: dict[str, str]) -> pd.DataFrame:
    table = pd.read_csv(META, sep="\t", dtype=str)
    table = table[["sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis"]].copy()
    table["donor"] = table["sample"].map(mapping).fillna(table["sample"])
    table["exclude"] = (
        table["exclude_stage_analysis"].str.lower().isin({"true", "t", "1", "yes"})
    )
    rows = []
    for donor, group in table.groupby("donor", sort=False):
        datasets = sorted(set(group["dataset"].dropna()) - {""})
        stages = sorted(set(group["disease_stage_coarse"].dropna()) - {""})
        rows.append({"donor": donor,
                     "dataset": datasets[0] if datasets else None,
                     "stage": stages[0] if stages else None,
                     "exclude": bool(group["exclude"].any())})
    return pd.DataFrame(rows)


def fit(frame: pd.DataFrame, ycol: str) -> dict:
    """Frozen primary model: score ~ stage_ordinal + factor(dataset)."""
    dummies = pd.get_dummies(frame["dataset"], drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(frame)),
                         frame["stage_ordinal"].to_numpy(float),
                         dummies.to_numpy(float)])
    y = frame[ycol].to_numpy(float)
    if np.linalg.matrix_rank(X) != X.shape[1]:
        return {"beta": np.nan, "hc3_pvalue": np.nan, "max_leverage": np.nan}
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    df = X.shape[0] - X.shape[1]
    H = X @ XtX_inv @ X.T
    h = np.clip(np.diag(H), 0.0, 1 - 1e-12)
    omega = (resid / (1.0 - h)) ** 2
    V = XtX_inv @ (X.T @ (omega[:, None] * X)) @ XtX_inv
    se_hc3 = float(np.sqrt(V[1, 1]))
    t = float(beta[1]) / se_hc3
    return {"beta": float(beta[1]),
            "hc3_pvalue": float(2 * stats.t.sf(abs(t), df)),
            "max_leverage": float(h.max())}


def main() -> None:
    scores = pd.read_csv(RES / "03_donor_scores_stored_vs_corrected.tsv", sep="\t")
    meta = donor_metadata(donor_map())
    rows = []
    for module in sorted(scores["module"].unique()):
        s = scores[scores["module"] == module]
        f = s.merge(meta, on="donor", how="inner", validate="one_to_one")
        f = f[(~f["exclude"]) & f["stage"].isin(PRIMARY_STAGES)].copy()
        f["stage_ordinal"] = f["stage"].map(PRIMARY_STAGES).astype(float)
        full_st = fit(f, "stored")
        full_co = fit(f, "corrected")
        pct = []
        for d in f["donor"]:
            g = f[f["donor"] != d]
            if g["dataset"].nunique() != f["dataset"].nunique():
                continue  # dropping the last donor of a dataset changes the design
            a, b = fit(g, "stored"), fit(g, "corrected")
            if np.isfinite(a["beta"]) and a["beta"] != 0:
                pct.append(100.0 * b["beta"] / a["beta"])
        pct = np.array(pct, dtype=float)
        rows.append({
            "module": module,
            "n_donors": int(len(f)),
            "full_pct_retained": 100.0 * full_co["beta"] / full_st["beta"],
            "loo_min_pct_retained": float(np.nanmin(pct)),
            "loo_max_pct_retained": float(np.nanmax(pct)),
            "loo_median_pct_retained": float(np.nanmedian(pct)),
            "loo_n_fits": int(np.isfinite(pct).sum()),
            "full_corrected_hc3_p": full_co["hc3_pvalue"],
            "max_leverage": full_st["max_leverage"],
        })
    out = pd.DataFrame(rows)
    out.to_csv(RES / "07_donor_jackknife.tsv", sep="\t", index=False)

    # which donor carries the highest leverage
    m8 = scores[scores["module"] == 8].merge(meta, on="donor", how="inner")
    m8 = m8[(~m8["exclude"]) & m8["stage"].isin(PRIMARY_STAGES)].copy()
    m8["stage_ordinal"] = m8["stage"].map(PRIMARY_STAGES).astype(float)
    dummies = pd.get_dummies(m8["dataset"], drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(m8)), m8["stage_ordinal"].to_numpy(float),
                         dummies.to_numpy(float)])
    h = np.diag(X @ np.linalg.inv(X.T @ X) @ X.T)
    m8 = m8.assign(leverage=h).sort_values("leverage", ascending=False)
    print("donors per dataset / stage:")
    print(m8.groupby(["dataset", "stage"]).size().to_string())
    print("\ntop leverage donors (module-8 design, identical for every module):")
    print(m8[["donor", "dataset", "stage", "leverage"]].head(8).to_string(index=False))

    print("\nleave-one-donor-out %-stage-beta-retained:")
    show = out[out["module"].isin(FOCUS)].copy()
    show["label"] = show["module"].map(FOCUS)
    cols = ["module", "label", "full_pct_retained", "loo_min_pct_retained",
            "loo_median_pct_retained", "loo_max_pct_retained", "loo_n_fits"]
    print(show[cols].to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print("\n[jackknife] wrote 07_donor_jackknife.tsv", flush=True)


if __name__ == "__main__":
    main()
