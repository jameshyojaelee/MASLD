#!/usr/bin/env python3
"""
34_vu_substate_within_hepatocytes.py — Within-hepatocyte F3b axis bimodality test.

Re-tests F3a / F3b axis bimodality using ONLY hepatocyte-dominant spots
(≥80% hepatocyte abundance). Same Bayesian-Information-Criterion + dip
Monte-Carlo tests as the original 32_vu_substate_f3a_f3b.py.

If F3a/F3b bimodality is a hepatocyte/stellate compositional artifact,
restricting to hep-dominant spots should ABOLISH the bimodality.

Inputs:
  RNA-seq/results/granular_staging/vu_deconvolved_substate.csv

Outputs:
  RNA-seq/results/granular_staging/vu_within_hepatocyte_bimodality.csv
    columns: signature, n_hep_dom_spots, BIC_diff, sep_sigma, dip_p,
             pass_real_bimodality

SLURM: --partition=cpu --qos=interactive --mem=32G --cpus=4 --time=4:00:00
       env: spatial
"""

import pathlib
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GRAN_OUT = PROJECT_ROOT / "RNA-seq/results/granular_staging"

DECONV_CSV = GRAN_OUT / "vu_deconvolved_substate.csv"
OUT_CSV = GRAN_OUT / "vu_within_hepatocyte_bimodality.csv"

# Pre-registered acceptance thresholds (from team M3 spec):
BIC_DIFF_THRESHOLD = 10.0
DIP_P_THRESHOLD = 0.05


def fit_gmm_2(x, n_iter=200, tol=1e-6):
    n = len(x)
    ord_ = np.sort(x)
    half = n // 2
    mu1, mu2 = ord_[:half].mean(), ord_[half:].mean()
    sd1, sd2 = max(ord_[:half].std(), 1e-3), max(ord_[half:].std(), 1e-3)
    p1 = 0.5
    ll_old = -np.inf
    for _ in range(n_iter):
        d1 = p1 * stats.norm.pdf(x, mu1, sd1)
        d2 = (1 - p1) * stats.norm.pdf(x, mu2, sd2)
        g1 = d1 / (d1 + d2 + 1e-300)
        p1 = g1.mean()
        mu1 = (g1 * x).sum() / g1.sum()
        mu2 = ((1 - g1) * x).sum() / (1 - g1).sum()
        sd1 = max(np.sqrt((g1 * (x - mu1) ** 2).sum() / g1.sum()), 1e-3)
        sd2 = max(np.sqrt(((1 - g1) * (x - mu2) ** 2).sum() / (1 - g1).sum()),
                  1e-3)
        ll = np.log(p1 * stats.norm.pdf(x, mu1, sd1) +
                    (1 - p1) * stats.norm.pdf(x, mu2, sd2) + 1e-300).sum()
        if abs(ll - ll_old) < tol:
            break
        ll_old = ll
    return dict(mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2, p1=p1, ll=ll)


def bic_compare(x):
    n = len(x)
    ll1 = stats.norm.logpdf(x, x.mean(), x.std()).sum()
    bic1 = -2 * ll1 + 2 * np.log(n)
    fit = fit_gmm_2(x)
    bic2 = -2 * fit["ll"] + 5 * np.log(n)
    sep = abs(fit["mu1"] - fit["mu2"]) / np.sqrt(
        (fit["sd1"] ** 2 + fit["sd2"] ** 2) / 2
    )
    return dict(bic1=bic1, bic2=bic2, bic_diff=bic1 - bic2,
                sep_sigma=sep, p1=fit["p1"])


def hartigan_mc(x, n_mc=999, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)
    grid = np.linspace(x.min(), x.max(), 200)
    e = np.array([np.mean(x <= g) for g in grid])
    g_pdf = stats.norm.cdf(grid, x.mean(), x.std())
    obs_dip = float(np.max(np.abs(e - g_pdf)))
    null_dips = []
    sx, mx = x.std(), x.mean()
    for _ in range(n_mc):
        y = rng.normal(mx, sx, len(x))
        ey = np.array([np.mean(y <= g) for g in grid])
        gy = stats.norm.cdf(grid, y.mean(), y.std())
        null_dips.append(np.max(np.abs(ey - gy)))
    p = (np.sum(np.array(null_dips) >= obs_dip) + 1) / (n_mc + 1)
    return obs_dip, p


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hep-threshold", type=float, default=0.80,
                        help="Hepatocyte fraction threshold")
    parser.add_argument("--n-mc", type=int, default=999,
                        help="Monte-Carlo iterations for dip test")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if not DECONV_CSV.exists():
        print(f"  ERROR: deconvolution output not found: {DECONV_CSV}")
        print(f"  Run 33_vu_cell2location_deconvolution.py first")
        sys.exit(1)

    if OUT_CSV.exists() and not args.force:
        print(f"  Output exists: {OUT_CSV} — pass --force to rerun")
        sys.exit(0)

    print(f"== Within-hepatocyte F3a/F3b bimodality ==")
    print(f"  Loading: {DECONV_CSV}")
    df = pd.read_csv(DECONV_CSV)
    print(f"  Total spots: {len(df)}")

    # Hepatocyte-dominant filter
    hep_dom = df["hep_fraction"] >= args.hep_threshold
    df_hep = df[hep_dom].copy()
    print(f"  Hep-dominant (≥{args.hep_threshold}): {len(df_hep)} "
          f"({100 * hep_dom.mean():.1f}%)")
    print(f"  Method used: {df['method'].iloc[0]}")

    # Re-z-score WITHIN hepatocyte-dominant subset, per individual
    # (so the bimodality test is comparable to original whole-tissue test)
    for sig in ["li_f3a", "li_f3b"]:
        df_hep[f"{sig}_z_hep"] = df_hep.groupby("individual")[sig].transform(
            lambda x: (x - x.mean()) / (x.std() + 1e-6)
        )
    df_hep["mixing_axis_hep"] = df_hep["li_f3b_z_hep"] - df_hep["li_f3a_z_hep"]

    # Run bimodality tests on each signature
    rng = np.random.default_rng(42)
    rows = []
    sigs = ["li_f3a", "li_f3b", "li_f3a_z_hep", "li_f3b_z_hep",
            "mixing_axis_hep", "mixing_axis"]
    for sig in sigs:
        if sig not in df_hep.columns:
            continue
        x = df_hep[sig].values.astype(float)
        x = x[~np.isnan(x)]
        if len(x) < 50:
            print(f"  {sig}: skipped (n={len(x)} too small)")
            continue
        bic = bic_compare(x)
        obs_dip, dip_p = hartigan_mc(x, n_mc=args.n_mc, rng=rng)
        passes = (bic["bic_diff"] >= BIC_DIFF_THRESHOLD) and (
            dip_p < DIP_P_THRESHOLD
        )
        rows.append({
            "signature": sig,
            "n_hep_dom_spots": len(x),
            "BIC_diff": bic["bic_diff"],
            "sep_sigma": bic["sep_sigma"],
            "GMM_p1": bic["p1"],
            "GMM_mu1": bic.get("mu1") if False else None,
            "dip_obs": obs_dip,
            "dip_p": dip_p,
            "pass_real_bimodality": passes,
        })
        flag = "PASS" if passes else "FAIL"
        print(f"  {sig:24s} n={len(x):6d} ΔBIC={bic['bic_diff']:9.2f} "
              f"sep={bic['sep_sigma']:5.2f}σ dip-p={dip_p:.4f} → {flag}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"  Saved: {OUT_CSV}")

    # Print summary
    n_pass = out["pass_real_bimodality"].sum()
    print(f"\n  Verdict: {n_pass}/{len(out)} signatures show real "
          f"within-hepatocyte bimodality "
          f"(ΔBIC≥{BIC_DIFF_THRESHOLD} AND dip-p<{DIP_P_THRESHOLD})")


if __name__ == "__main__":
    main()
