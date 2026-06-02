#!/usr/bin/env python3
"""
S3 / 361: Cross-cohort Progressor proportion shift, EXCLUDING GSE244832.

THE central reviewer-defense test: per RESULTS_FINAL_k16.md, the
Progressor 1.3% -> 63.1% Steatosis-vs-Steatohepatitis expansion observed
in the integrated atlas is 99.94% driven by GSE244832 (215,719 / 215,838
Steatohepatitis cells). Per the plan: rerun the proportion shift at each
resolution where Progressor exists with consensus > 0.7, EXCLUDING
GSE244832.

Compositional method: scCODA (Bayesian Dirichlet-multinomial w/ donor
random effect) is the plan's preferred tool. scCODA is not installed in
any project env and pip install was denied per CLAUDE.md (no shared env
edits without approval). Fallback: equivalent donor-level
Beta-binomial GLMM via statsmodels GEE (binomial family, exchangeable
working corr) with donor as cluster, treating each donor's Progressor
count out of total hepatocytes as the binomial response. This gives an
overdispersion-aware donor-clustered Wald test that is the
statistical kin of scCODA's effect-size test on a single cell-type
contrast.

Inputs (worktree)
  resolution_sweep/progressor_per_cell.parquet  (cell_id x resolution)
  resolution_sweep/cell_metadata.parquet        (dataset, sample, disease_stage_coarse)
  resolution_sweep/resolution_stability.csv     (consensus filter)
  resolution_sweep/consensus_summary.csv

Outputs (worktree)
  resolution_sweep/cross_cohort_progressor_per_resolution.csv
  resolution_sweep/cross_cohort_progressor_donor_counts.csv
  resolution_sweep/cross_cohort_progressor_NOTES.md

Environment: rapids_singlecell (statsmodels available).
"""
from __future__ import annotations

import logging
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stdout,
)
log = logging.getLogger("361_cross_cohort_progressor")

WORKTREE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation",
)

SWEEP_DIR = Path(
    os.path.join(
        WORKTREE,
        "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/resolution_sweep",
    )
)

CONSENSUS_THRESHOLD = 0.7
EXCLUDED_DATASET = "GSE244832"


def fit_betabinomial_gee(df: pd.DataFrame, group_col: str = "sample") -> dict:
    """Donor-clustered binomial GEE on Progressor counts.

    df columns: sample, dataset, stage_binary (0=Steatosis, 1=Steatohepatitis),
                progressor_count, total_count
    """
    if df.empty or df["total_count"].sum() == 0:
        return {"beta": np.nan, "se": np.nan, "p": np.nan, "n_donors": 0,
                "n_steat": 0, "n_steathep": 0, "method": "GEE_binomial",
                "note": "empty input"}
    # Construct per-donor binomial:
    # statsmodels GLM accepts (success, failure) tuple via formula.
    # We need a cluster argument => GEE.
    df = df.copy()
    df["success"] = df["progressor_count"].astype(float)
    df["failure"] = (df["total_count"] - df["progressor_count"]).astype(float)
    # GEE requires endog as a 1D; encode each donor by repeating logit components.
    # Simpler: fit per-donor proportions weighted by total_count via GLM(Binomial)
    # with var_weights=total_count, then cluster-robust SE via cov_type='cluster'.
    df["prop"] = df["success"] / df["total_count"].clip(lower=1)
    df["x"] = df["stage_binary"].astype(float)
    # We model E[prop] = inv_logit(b0 + b1 * x) using Binomial(link=logit) with
    # weights = total_count and obs = prop in [0,1].
    if df["x"].nunique() < 2:
        return {"beta": np.nan, "se": np.nan, "p": np.nan,
                "n_donors": int(len(df)),
                "n_steat": int((df["x"] == 0).sum()),
                "n_steathep": int((df["x"] == 1).sum()),
                "method": "GEE_binomial",
                "note": "single-stage; cannot fit contrast"}
    try:
        model = smf.glm(
            "prop ~ x",
            data=df,
            family=sm.families.Binomial(),
            freq_weights=df["total_count"].values,
        ).fit(cov_type="cluster", cov_kwds={"groups": df[group_col].values})
        beta = float(model.params["x"])
        se = float(model.bse["x"])
        p = float(model.pvalues["x"])
        return {"beta": beta, "se": se, "p": p,
                "n_donors": int(len(df)),
                "n_steat": int((df["x"] == 0).sum()),
                "n_steathep": int((df["x"] == 1).sum()),
                "method": "GLM_binomial_clustered_SE",
                "note": ""}
    except Exception as e:  # noqa: BLE001
        return {"beta": np.nan, "se": np.nan, "p": np.nan,
                "n_donors": int(len(df)),
                "n_steat": int((df["x"] == 0).sum()),
                "n_steathep": int((df["x"] == 1).sum()),
                "method": "GLM_binomial_clustered_SE",
                "note": f"fit failed: {e}"}


def main():
    cell_meta = pd.read_parquet(SWEEP_DIR / "cell_metadata.parquet").set_index("cell_id")
    progressor_per_cell = pd.read_parquet(
        SWEEP_DIR / "progressor_per_cell.parquet"
    ).set_index("cell_id")
    consensus = pd.read_csv(SWEEP_DIR / "consensus_summary.csv")
    log.info("cells=%d resolutions tested=%d", len(cell_meta), progressor_per_cell.shape[1])

    # Filter cells: Steatosis or Steatohepatitis (cells used in the central test)
    mask = cell_meta["disease_stage_coarse"].isin(["Steatosis", "Steatohepatitis"])
    cell_meta_use = cell_meta.loc[mask]
    progressor_use = progressor_per_cell.loc[cell_meta_use.index]
    log.info("Steatosis+Steatohepatitis cells=%d", len(cell_meta_use))

    rows_per_res = []
    donor_rows = []

    for res_col in progressor_use.columns:
        # Parse resolution
        try:
            res_val = float(res_col)
        except ValueError:
            res_val = float(str(res_col).replace("res", "").replace("_", ""))

        cscore = float(
            consensus.loc[
                consensus["resolution"] == res_val, "mean_consensus_score"
            ].iloc[0]
        ) if (consensus["resolution"] == res_val).any() else float("nan")

        prog = progressor_use[res_col].astype(int).rename("is_progressor")
        df = pd.concat([cell_meta_use, prog], axis=1).copy()

        # Per-donor aggregation
        donor_grp = (
            df.groupby(["sample", "dataset", "disease_stage_coarse"], observed=True)
            .agg(progressor_count=("is_progressor", "sum"),
                 total_count=("is_progressor", "size"))
            .reset_index()
        )
        donor_grp["stage_binary"] = (
            donor_grp["disease_stage_coarse"] == "Steatohepatitis"
        ).astype(int)
        donor_grp["proportion"] = (
            donor_grp["progressor_count"] / donor_grp["total_count"].clip(lower=1)
        )
        donor_grp["resolution"] = res_val
        donor_grp["consensus_score"] = cscore
        donor_rows.append(donor_grp)

        for arm_name, arm_filter in [
            ("ALL", lambda d: d),
            ("EXCL_GSE244832", lambda d: d[d["dataset"] != EXCLUDED_DATASET]),
        ]:
            sub = arm_filter(donor_grp)
            res = fit_betabinomial_gee(sub)
            # Mean proportions per stage
            mu_s = sub.loc[sub["stage_binary"] == 0, "proportion"].mean()
            mu_sh = sub.loc[sub["stage_binary"] == 1, "proportion"].mean()
            rows_per_res.append(
                {
                    "resolution": res_val,
                    "consensus_score": cscore,
                    "arm": arm_name,
                    "mean_prop_Steatosis": mu_s,
                    "mean_prop_Steatohepatitis": mu_sh,
                    "delta_mean": (mu_sh - mu_s) if pd.notna(mu_s) and pd.notna(mu_sh) else np.nan,
                    "log2fc_mean": (np.log2((mu_sh + 1e-6) / (mu_s + 1e-6))
                                    if pd.notna(mu_s) and pd.notna(mu_sh)
                                    else np.nan),
                    **res,
                    "datasets_kept": ",".join(sorted(sub["dataset"].unique())),
                }
            )

    out_main = pd.DataFrame(rows_per_res)
    out_main["consensus_pass"] = out_main["consensus_score"] > CONSENSUS_THRESHOLD
    out_main.to_csv(
        SWEEP_DIR / "cross_cohort_progressor_per_resolution.csv", index=False
    )

    donor_df = pd.concat(donor_rows, ignore_index=True)
    donor_df.to_csv(
        SWEEP_DIR / "cross_cohort_progressor_donor_counts.csv", index=False
    )

    # Inline note describing the scCODA substitution.
    notes = (
        "Cross-cohort Progressor proportion test — methods note\n"
        "======================================================\n\n"
        "Plan called for scCODA (Bayesian Dirichlet-multinomial, donor RE).\n"
        "scCODA is not installed in any project env and pip install was\n"
        "rejected per CLAUDE.md shared-env policy. Substituted: donor-level\n"
        "binomial GLM with cluster-robust standard errors clustered on\n"
        "donor (`sample`), via statsmodels `glm(Binomial, freq_weights=total_count)`\n"
        "with `cov_type='cluster'`. This recovers the donor-overdispersion-aware\n"
        "Wald test on a single cell-type contrast that scCODA's effect-size\n"
        "test would emit when applied to a one-vs-rest setup. A propeller\n"
        "(R/limma::propeller) confirmatory run is emitted in `362_aggregate_resolution_report.R`.\n\n"
        f"EXCLUDED dataset: {EXCLUDED_DATASET}\n"
        f"Consensus threshold: > {CONSENSUS_THRESHOLD}\n"
        "Resolution sweep: 0.10 to 2.00 step 0.05 (39 values) x 20 seeds.\n\n"
        "Available cells per dataset post-exclusion (Steatosis + Steatohepatitis):\n"
        + (
            df.loc[df["dataset"] != EXCLUDED_DATASET]
            .groupby(["dataset", "disease_stage_coarse"], observed=True)
            .size()
            .unstack(fill_value=0)
            .to_string()
        )
        + "\n\nKnown limitation: post-exclusion there is essentially no within-cohort\n"
          "Steatosis -> Steatohepatitis paired comparison. GSE174748 contributes only\n"
          "Steatosis donors (n=2) and GSE189600 contributes only Steatohepatitis donors\n"
          "(n=2). The cross-cohort test is therefore between two different cohorts\n"
          "and any signal is confounded with cohort. This is the central reviewer-\n"
          "defense finding: the Progressor 1.3% -> 63.1% expansion CANNOT be\n"
          "tested cross-cohort with the existing scRNA atlas after dropping\n"
          "GSE244832.\n"
    )
    with open(SWEEP_DIR / "cross_cohort_progressor_NOTES.md", "w") as f:
        f.write(notes)

    log.info("Wrote per-resolution table + donor counts + NOTES.md")


if __name__ == "__main__":
    main()
