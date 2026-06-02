#!/usr/bin/env python3
"""
35_vu_moran_matched_null.py — Matched-null Moran's I for F3a/F3b spatial structure.

Tests whether observed Moran's I for the F3a / F3b / mixing-axis signatures
is in the >=80th percentile of a null distribution of 100 random gene sets
matched on size + mean expression + variance. Computed per individual
(n=5 patients), not pooled.

Also runs a per-patient mixed-effects (random-intercept) model on Moran's I
to assess between-patient consistency.

Inputs:
  - Vu spatial AnnData: Analysis/Spatial/results/preprocessed/merged_spatial_vu.h5ad
  - F3a / F3b signatures (Li et al., from script 32)

Outputs:
  - RNA-seq/results/granular_staging/vu_moran_matched_null.csv
    columns: patient, signature, morans_I_observed, null_mean, null_sd,
             null_pct, morans_I_null_p, pass_real_organization
  - RNA-seq/results/granular_staging/vu_per_patient_mixed_effects.csv
    columns: signature, n_patients, fixed_effect_intercept, fe_se, fe_z,
             fe_pvalue, var_intercept, ci_lo, ci_hi, ci_excludes_zero,
             n_patients_individual_ci_excludes_zero

SLURM: --partition=cpu --qos=interactive --mem=64G --cpus=8 --time=12:00:00
       env: spatial
"""

import pathlib
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SPATIAL_DIR = PROJECT_ROOT / "Analysis/Spatial"
GRAN_OUT = PROJECT_ROOT / "RNA-seq/results/granular_staging"

VU_NORM_H5AD = SPATIAL_DIR / "results/preprocessed/merged_spatial_vu.h5ad"
OUT_MORAN_CSV = GRAN_OUT / "vu_moran_matched_null.csv"
OUT_MIXED_CSV = GRAN_OUT / "vu_per_patient_mixed_effects.csv"

LI_F3A = [
    "HSPA5", "DDIT3", "ATF4", "ATF6", "EIF2AK3", "ERN1", "XBP1", "DNAJB9",
    "SREBF1", "SREBF2", "HMGCR", "HMGCS1", "SQLE", "DHCR7", "DHCR24",
    "MVD", "MVK", "FASN", "ACACA", "SCD", "INSIG1",
    "FABP1", "ACOX1", "CPT1A", "PLIN2", "DGAT1",
]
LI_F3B = [
    "IGFBP7", "BGN", "COL1A2", "COL3A1", "TIMP1",
    "COL1A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "FN1", "VCAN", "DCN", "LUM", "LOX", "LOXL1", "LOXL2",
    "MMP2", "MMP9", "TIMP2", "ACTA2", "TAGLN",
    "CDKN1A", "CDKN2A", "GLB1", "SERPINE1",
]


def get_spot_score(adata_sub, gene_list, rng=None, score_name="_score"):
    """Score sub-AnnData on gene_list, return per-spot score vector."""
    present = [g for g in gene_list if g in adata_sub.var_names]
    if len(present) < 3:
        return None
    sc.tl.score_genes(adata_sub, gene_list=present, score_name=score_name,
                      random_state=42)
    return adata_sub.obs[score_name].values.astype(float)


def matched_random_genes(adata, target_genes, n_genes, rng):
    """Return random gene set matched to target on size + mean + variance.

    Matched within mean-expression decile and variance decile bins.
    """
    target_present = [g for g in target_genes if g in adata.var_names]
    if not target_present:
        return None
    # Pre-compute gene mean / var if not cached
    if "_gene_mean" not in adata.var.columns:
        X = adata.X
        if hasattr(X, "toarray"):
            X_arr = X.toarray()
        else:
            X_arr = np.asarray(X)
        adata.var["_gene_mean"] = X_arr.mean(axis=0)
        adata.var["_gene_var"] = X_arr.var(axis=0)

    target_mean = adata.var.loc[target_present, "_gene_mean"].mean()
    target_var = adata.var.loc[target_present, "_gene_var"].mean()
    # Define eligibility window: ±50% mean, ±50% var
    mean_lo, mean_hi = target_mean * 0.5, target_mean * 1.5
    var_lo, var_hi = target_var * 0.5, target_var * 1.5
    eligible = adata.var[
        (adata.var["_gene_mean"] >= mean_lo)
        & (adata.var["_gene_mean"] <= mean_hi)
        & (adata.var["_gene_var"] >= var_lo)
        & (adata.var["_gene_var"] <= var_hi)
    ].index.tolist()
    eligible = [g for g in eligible if g not in set(target_genes)]

    if len(eligible) < n_genes:
        # Relax: drop variance constraint, keep mean only
        eligible = adata.var[
            (adata.var["_gene_mean"] >= mean_lo)
            & (adata.var["_gene_mean"] <= mean_hi)
        ].index.tolist()
        eligible = [g for g in eligible if g not in set(target_genes)]
    if len(eligible) < n_genes:
        eligible = [g for g in adata.var_names if g not in set(target_genes)]
    return list(rng.choice(eligible, size=n_genes, replace=False))


def compute_morans_i(values, coords, k=6):
    """Per-spot Moran's I given values & coordinates."""
    from libpysal.weights import KNN
    from esda.moran import Moran
    valid = ~np.isnan(values)
    if valid.sum() < 50:
        return np.nan, np.nan
    kw = KNN.from_array(coords[valid], k=k)
    mor = Moran(values[valid], kw, permutations=99)
    return mor.I, mor.p_sim


def per_patient_morans_with_null(adata, sig_name, target_genes,
                                  n_null=100, rng=None):
    """For each patient: compute observed Moran's I + matched-null distribution.

    Returns list of dicts (one per patient).
    """
    if rng is None:
        rng = np.random.default_rng(42)
    rows = []
    patients = sorted(adata.obs["individual"].unique())
    for ind in patients:
        sub = adata[adata.obs["individual"] == ind].copy()
        if sub.n_obs < 100 or "spatial" not in sub.obsm:
            print(f"    {ind}: skipped (n={sub.n_obs})")
            continue
        coords = sub.obsm["spatial"]

        # Observed
        obs_score = get_spot_score(sub, target_genes,
                                   score_name=f"_{sig_name}_obs")
        if obs_score is None:
            print(f"    {ind}: signature {sig_name} unavailable")
            continue
        obs_I, obs_p = compute_morans_i(obs_score, coords, k=6)

        # Null distribution
        null_Is = []
        for nn in range(n_null):
            null_genes = matched_random_genes(
                adata, target_genes, n_genes=len(target_genes), rng=rng
            )
            null_score = get_spot_score(
                sub, null_genes,
                score_name=f"_{sig_name}_null_{nn}",
            )
            if null_score is None:
                continue
            null_I, _ = compute_morans_i(null_score, coords, k=6)
            null_Is.append(null_I)
            # Drop temp obs col
            try:
                sub.obs.drop(columns=[f"_{sig_name}_null_{nn}"], inplace=True)
            except KeyError:
                pass

        null_Is = np.array([x for x in null_Is if not np.isnan(x)])
        if len(null_Is) < 10:
            print(f"    {ind}/{sig_name}: too few null (n={len(null_Is)})")
            continue
        null_mean = null_Is.mean()
        null_sd = null_Is.std()
        # Empirical p (right-tailed)
        emp_p = (np.sum(null_Is >= obs_I) + 1) / (len(null_Is) + 1)
        # Percentile of obs
        null_pct = (np.sum(null_Is < obs_I) / len(null_Is)) * 100.0
        passes = null_pct >= 80.0
        rows.append({
            "patient": ind,
            "signature": sig_name,
            "morans_I_observed": obs_I,
            "morans_I_obs_perm_p": obs_p,
            "n_null": len(null_Is),
            "null_mean": null_mean,
            "null_sd": null_sd,
            "null_pct": null_pct,
            "morans_I_null_p": emp_p,
            "pass_real_organization": passes,
        })
        flag = "PASS" if passes else "FAIL"
        print(f"    {ind:8s} {sig_name:18s} obs={obs_I:6.3f} "
              f"null={null_mean:6.3f}±{null_sd:5.3f} "
              f"pct={null_pct:5.1f} → {flag}")
    return rows


def fit_mixed_effects(per_patient_df, sig_name):
    """Random-intercept mixed-effects model on Moran's I per patient.

    Tests whether the patient-averaged Moran's I differs from 0 with a
    random intercept per patient.

    With only ~5 patients, this is essentially a one-sample t-test +
    random-effect variance. Use a one-sample t-test against 0 as the
    primary test, since 5 random levels is too small to reliably estimate
    a variance component.
    """
    sub = per_patient_df[per_patient_df["signature"] == sig_name].copy()
    if len(sub) < 2:
        return None

    vals = sub["morans_I_observed"].values.astype(float)
    n = len(vals)
    mean = vals.mean()
    sd = vals.std(ddof=1) if n > 1 else 0.0
    se = sd / np.sqrt(n) if n > 1 else 0.0
    if se > 0:
        from scipy.stats import t as stt
        tstat = mean / se
        pval = 2 * (1 - stt.cdf(abs(tstat), df=n - 1))
        # 95% CI
        tcrit = stt.ppf(0.975, df=n - 1)
        ci_lo = mean - tcrit * se
        ci_hi = mean + tcrit * se
    else:
        tstat, pval, ci_lo, ci_hi = np.nan, np.nan, np.nan, np.nan

    # Per-patient permutation p<0.05 count
    n_pass_individual = (sub["morans_I_obs_perm_p"] < 0.05).sum()

    return {
        "signature": sig_name,
        "n_patients": n,
        "fixed_effect_intercept": mean,
        "fe_se": se,
        "fe_t": tstat,
        "fe_pvalue": pval,
        "var_intercept": float(np.var(vals, ddof=1)) if n > 1 else 0.0,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "ci_excludes_zero": (ci_lo > 0) or (ci_hi < 0),
        "n_patients_individual_perm_p_lt_05": int(n_pass_individual),
        "pass_4_of_5_perm": int(n_pass_individual >= 4),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-null", type=int, default=100,
                        help="Number of matched-null gene sets per signature")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if OUT_MORAN_CSV.exists() and not args.force:
        print(f"  Output exists: {OUT_MORAN_CSV} — pass --force to rerun")
        sys.exit(0)

    print(f"== Matched-null Moran's I (per patient) ==")
    print(f"  Loading: {VU_NORM_H5AD}")
    adata = sc.read_h5ad(VU_NORM_H5AD)
    print(f"  Spots: {adata.n_obs}, Genes: {adata.n_vars}")
    patients = sorted(adata.obs["individual"].unique())
    print(f"  Patients: {len(patients)} ({patients})")

    signatures = {
        "li_f3a": LI_F3A,
        "li_f3b": LI_F3B,
        # mixing axis: F3b - F3a; signed combination of two gene sets,
        # treat as their union for null-matching purposes
        "li_f3b_minus_f3a": ("MIX", LI_F3A, LI_F3B),
    }

    rng = np.random.default_rng(42)
    all_rows = []

    for sig_name, payload in signatures.items():
        print(f"\n  Signature: {sig_name}")
        if isinstance(payload, tuple) and payload[0] == "MIX":
            f3a_genes, f3b_genes = payload[1], payload[2]
            # For mixing axis: score F3a + F3b separately, take F3b_z - F3a_z
            # PER PATIENT (since z is per-patient).
            rows = []
            for ind in patients:
                sub = adata[adata.obs["individual"] == ind].copy()
                if sub.n_obs < 100 or "spatial" not in sub.obsm:
                    continue
                coords = sub.obsm["spatial"]
                f3a_score = get_spot_score(sub, f3a_genes, score_name="_mix_f3a")
                f3b_score = get_spot_score(sub, f3b_genes, score_name="_mix_f3b")
                if f3a_score is None or f3b_score is None:
                    continue
                f3a_z = (f3a_score - f3a_score.mean()) / (f3a_score.std() + 1e-6)
                f3b_z = (f3b_score - f3b_score.mean()) / (f3b_score.std() + 1e-6)
                obs_score = f3b_z - f3a_z
                obs_I, obs_p = compute_morans_i(obs_score, coords, k=6)

                null_Is = []
                target_combined = list(set(f3a_genes) | set(f3b_genes))
                for nn in range(args.n_null):
                    null_a = matched_random_genes(adata, f3a_genes,
                                                  len(f3a_genes), rng)
                    null_b = matched_random_genes(adata, f3b_genes,
                                                  len(f3b_genes), rng)
                    null_a_score = get_spot_score(sub, null_a,
                                                  score_name=f"_null_a_{nn}")
                    null_b_score = get_spot_score(sub, null_b,
                                                  score_name=f"_null_b_{nn}")
                    if null_a_score is None or null_b_score is None:
                        continue
                    null_a_z = ((null_a_score - null_a_score.mean())
                                / (null_a_score.std() + 1e-6))
                    null_b_z = ((null_b_score - null_b_score.mean())
                                / (null_b_score.std() + 1e-6))
                    null_score = null_b_z - null_a_z
                    null_I, _ = compute_morans_i(null_score, coords, k=6)
                    null_Is.append(null_I)
                    try:
                        sub.obs.drop(columns=[f"_null_a_{nn}",
                                              f"_null_b_{nn}"],
                                     inplace=True, errors="ignore")
                    except (KeyError, TypeError):
                        pass

                null_Is = np.array([x for x in null_Is if not np.isnan(x)])
                if len(null_Is) < 10:
                    continue
                null_mean = null_Is.mean()
                null_sd = null_Is.std()
                # Two-tailed for mixing axis (sign matters)
                emp_p_right = (np.sum(null_Is >= obs_I) + 1) / (len(null_Is) + 1)
                emp_p_left = (np.sum(null_Is <= obs_I) + 1) / (len(null_Is) + 1)
                emp_p = 2 * min(emp_p_right, emp_p_left)
                null_pct = (np.sum(null_Is < obs_I) / len(null_Is)) * 100.0
                passes = (null_pct >= 80.0) or (null_pct <= 20.0)
                rows.append({
                    "patient": ind,
                    "signature": sig_name,
                    "morans_I_observed": obs_I,
                    "morans_I_obs_perm_p": obs_p,
                    "n_null": len(null_Is),
                    "null_mean": null_mean,
                    "null_sd": null_sd,
                    "null_pct": null_pct,
                    "morans_I_null_p": emp_p,
                    "pass_real_organization": passes,
                })
                flag = "PASS" if passes else "FAIL"
                print(f"    {ind:8s} {sig_name:18s} obs={obs_I:6.3f} "
                      f"null={null_mean:6.3f}±{null_sd:5.3f} "
                      f"pct={null_pct:5.1f} → {flag}")
            all_rows.extend(rows)
        else:
            target_genes = payload
            rows = per_patient_morans_with_null(
                adata, sig_name, target_genes,
                n_null=args.n_null, rng=rng,
            )
            all_rows.extend(rows)

    moran_df = pd.DataFrame(all_rows)
    moran_df.to_csv(OUT_MORAN_CSV, index=False)
    print(f"\n  Saved: {OUT_MORAN_CSV} ({len(moran_df)} rows)")

    # Mixed-effects per signature
    print(f"\n== Per-patient mixed-effects (random intercept) ==")
    me_rows = []
    for sig in moran_df["signature"].unique():
        me = fit_mixed_effects(moran_df, sig)
        if me is None:
            continue
        me_rows.append(me)
        print(f"  {sig:20s} n={me['n_patients']} "
              f"FE={me['fixed_effect_intercept']:6.3f} "
              f"95% CI=[{me['ci_lo']:6.3f}, {me['ci_hi']:6.3f}] "
              f"p={me['fe_pvalue']:.4f} "
              f"individual_p<.05: {me['n_patients_individual_perm_p_lt_05']}/{me['n_patients']}")

    me_df = pd.DataFrame(me_rows)
    me_df.to_csv(OUT_MIXED_CSV, index=False)
    print(f"\n  Saved: {OUT_MIXED_CSV}")


if __name__ == "__main__":
    main()
