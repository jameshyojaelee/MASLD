#!/usr/bin/env python3
"""
153_pseudo_survival.py — Pseudo-Survival on Pseudotime
-------------------------------------------------------
Cox-like proportional hazards analysis using molecular pseudotime as
"time" and advanced fibrosis (F>=3) as "event". Since lifelines and
scikit-survival are not installed, we implement:

1. **Logistic regression** of event ~ features (odds ratios as hazard
   ratio proxies).
2. **Pseudo-C-index**: concordance between predicted risk and pseudotime
   ordering among event-positive samples, approximating survival
   concordance.
3. **Stratified survival curves**: CPS-quartile-stratified "survival"
   (fraction without advanced fibrosis) across pseudotime bins.

CRITICAL CAVEAT (must be stated in any manuscript use):
  Pseudotime is NOT calendar time. This is molecular trajectory analysis,
  not true survival prediction. All results should be framed as:
  "features associated with advanced disease position in molecular
  trajectory space." No temporal precedence can be established from
  cross-sectional data. True survival analysis requires longitudinal
  cohorts.

Input files:
  - results/prognosis/prognosis_feature_matrix.csv
  - results/prognosis/prognosis_pseudo_labels.csv
  - results/progression/consensus_pseudotime.csv

Output (to results/prognosis/):
  - survival_hazard_ratios.csv        (feature odds ratios + p-values)
  - survival_curves_by_cps.csv        (CPS-quartile survival curves)
  - pseudo_survival_c_index.csv       (concordance metrics)

SLURM: io, 4 CPUs, 8G RAM, 48h
Env:   micromamba activate spatial
"""

import logging
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")

FEATURE_MATRIX = os.path.join(INTEG, "results/prognosis/prognosis_feature_matrix.csv")
PSEUDO_LABELS = os.path.join(INTEG, "results/prognosis/prognosis_pseudo_labels.csv")
CONSENSUS_PT = os.path.join(INTEG, "results/progression/consensus_pseudotime.csv")

OUT_DIR = os.path.join(INTEG, "results/prognosis")

# ── Parameters ───────────────────────────────────────────────────────────────
N_PSEUDOTIME_BINS = 20  # bins for survival curves
CPS_QUARTILE_LABELS = ["Q1_low", "Q2", "Q3", "Q4_high"]
SEED = 42


# ── Utility ──────────────────────────────────────────────────────────────────
def safe_load(path, **kwargs):
    """Load CSV with existence check."""
    if not os.path.exists(path):
        log.warning("File not found: %s", path)
        return None
    log.info("Loading %s", os.path.basename(path))
    return pd.read_csv(path, **kwargs)


def concordance_index(event_times, predicted_risk, event_observed):
    """
    Compute Harrell's C-index (concordance index).

    For each concordant pair (i, j) where event_time_i < event_time_j and
    event_i == 1: check if predicted_risk_i > predicted_risk_j.

    Returns: c_index, n_concordant, n_discordant, n_tied, n_pairs
    """
    n = len(event_times)
    concordant = 0
    discordant = 0
    tied = 0
    n_pairs = 0

    for i in range(n):
        if not event_observed[i]:
            continue  # only count pairs where the earlier sample has an event
        for j in range(n):
            if i == j:
                continue
            if event_times[j] <= event_times[i]:
                continue  # j must have later or equal time
            # Valid pair: i had event, j has later time
            n_pairs += 1
            if predicted_risk[i] > predicted_risk[j]:
                concordant += 1
            elif predicted_risk[i] < predicted_risk[j]:
                discordant += 1
            else:
                tied += 1

    if n_pairs == 0:
        return 0.5, 0, 0, 0, 0  # no valid pairs

    c_index = (concordant + 0.5 * tied) / n_pairs
    return c_index, concordant, discordant, tied, n_pairs


def concordance_index_fast(event_times, predicted_risk, event_observed):
    """
    Vectorized approximate C-index for large datasets.
    Uses a sampling approach if n_pairs would be very large.
    """
    n = len(event_times)
    # For datasets of this size (~1000 samples), direct computation is fine
    if n > 2000:
        # Subsample for speed
        rng = np.random.RandomState(SEED)
        idx = rng.choice(n, min(n, 2000), replace=False)
        event_times = event_times[idx]
        predicted_risk = predicted_risk[idx]
        event_observed = event_observed[idx]

    return concordance_index(event_times, predicted_risk, event_observed)


# ══════════════════════════════════════════════════════════════════════════════
# PART 1: Feature hazard ratios (logistic regression approximation)
# ══════════════════════════════════════════════════════════════════════════════
def compute_hazard_ratios(features, labels, pseudotime):
    """
    Logistic regression of F>=3 event ~ features.
    Reports odds ratios as hazard ratio proxies.
    """
    log.info("=" * 70)
    log.info("PART 1: Feature Hazard Ratios (Logistic Regression Proxy)")
    log.info("=" * 70)
    log.info("NOTE: These are odds ratios from logistic regression, reported as")
    log.info("hazard ratio proxies. True HRs require Cox PH (lifelines not installed).")
    log.info("")

    # Build analysis dataframe
    fib = pd.to_numeric(labels["fibrosis_stage"], errors="coerce")
    event = (fib >= 3).astype(float)
    event[fib.isna() | (fib < 0)] = np.nan

    # Select features: focused set
    feature_candidates = [
        "ct_Stellate", "ct_stellate_gt4pct", "ct_stellate_hepatocyte_ratio",
        "ct_fibrogenic_axis", "ct_immune_infiltration",
        "clin_sex",
        "trans_F0_to_F1", "trans_F1_to_F2", "trans_F2_to_F3", "trans_F3_to_F4",
        "hcc_genetic_score",
        "gas6_mertk_product", "gas6_mertk_ratio",
    ]
    # Add top 10 divergence genes
    div_cols = [c for c in features.columns if c.startswith("div_")][:10]
    feature_candidates.extend(div_cols)

    feature_cols = [c for c in feature_candidates if c in features.columns]
    log.info("Feature set: %d features", len(feature_cols))

    # Valid samples: non-NaN event, non-NaN pseudotime, non-NaN for most features
    valid = event.notna() & pseudotime.notna()
    n_valid = valid.sum()
    log.info("Samples with valid event + pseudotime: %d", n_valid)

    if n_valid < 50:
        log.warning("Too few valid samples (%d < 50). Skipping.", n_valid)
        return pd.DataFrame()

    X = features.loc[valid, feature_cols].values.astype(float)
    y = event[valid].values
    pt = pseudotime[valid].values

    # Impute NaN with column median
    for j in range(X.shape[1]):
        nan_mask = np.isnan(X[:, j])
        if nan_mask.any():
            med = np.nanmedian(X[:, j])
            X[nan_mask, j] = med if not np.isnan(med) else 0.0

    log.info("Event rate: %d / %d (%.1f%%)", y.sum(), len(y), 100 * y.mean())

    # --- Univariate odds ratios ---
    log.info("Computing univariate odds ratios ...")
    hr_rows = []
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    for j, feat in enumerate(feature_cols):
        x_j = X_scaled[:, j].reshape(-1, 1)
        lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=5000)
        lr.fit(x_j, y)
        coef = lr.coef_[0, 0]
        or_val = np.exp(coef)  # odds ratio per 1-SD increase

        # Wald test for significance
        proba = lr.predict_proba(x_j)[:, 1]
        W = proba * (1 - proba)
        # Fisher information
        fisher_info = np.sum(W * x_j[:, 0] ** 2)
        if fisher_info > 0:
            se = 1.0 / np.sqrt(fisher_info)
            z_stat = coef / se
            p_val = 2 * stats.norm.sf(np.abs(z_stat))
        else:
            se = np.nan
            z_stat = np.nan
            p_val = np.nan

        hr_rows.append({
            "feature": feat,
            "log_odds_ratio": coef,
            "odds_ratio": or_val,
            "se": se,
            "z_statistic": z_stat,
            "p_value": p_val,
            "or_lower_95": np.exp(coef - 1.96 * se) if not np.isnan(se) else np.nan,
            "or_upper_95": np.exp(coef + 1.96 * se) if not np.isnan(se) else np.nan,
            "model_type": "univariate_logistic",
            "note": "Odds ratio per 1-SD increase (standardized). Proxy for hazard ratio.",
        })

    # --- Multivariate model ---
    log.info("Fitting multivariate logistic regression ...")
    lr_multi = LogisticRegression(
        penalty="l2", C=1.0, solver="lbfgs", max_iter=5000, random_state=SEED
    )
    lr_multi.fit(X_scaled, y)
    multi_coefs = lr_multi.coef_[0]
    for j, feat in enumerate(feature_cols):
        coef = multi_coefs[j]
        or_val = np.exp(coef)
        hr_rows.append({
            "feature": feat,
            "log_odds_ratio": coef,
            "odds_ratio": or_val,
            "se": np.nan,  # Wald SE not trivial for regularized model
            "z_statistic": np.nan,
            "p_value": np.nan,
            "or_lower_95": np.nan,
            "or_upper_95": np.nan,
            "model_type": "multivariate_logistic_L2",
            "note": "Multivariate L2-regularized. Adjusted for all features.",
        })

    hr_df = pd.DataFrame(hr_rows)

    # Log top univariate results
    uni_df = hr_df[hr_df["model_type"] == "univariate_logistic"].copy()
    uni_df = uni_df.sort_values("p_value")
    log.info("Top 10 univariate features (by p-value):")
    for _, row in uni_df.head(10).iterrows():
        log.info("  %s: OR=%.3f [%.3f-%.3f], p=%.2e",
                 row["feature"], row["odds_ratio"],
                 row["or_lower_95"], row["or_upper_95"], row["p_value"])

    return hr_df


# ══════════════════════════════════════════════════════════════════════════════
# PART 2: Survival curves stratified by CPS quartile
# ══════════════════════════════════════════════════════════════════════════════
def survival_curves_by_cps(labels, pseudotime):
    """
    Kaplan-Meier-like curves: fraction of samples WITHOUT advanced fibrosis
    (F>=3) at each pseudotime bin, stratified by CPS quartile.
    """
    log.info("=" * 70)
    log.info("PART 2: Survival Curves by CPS Quartile")
    log.info("=" * 70)

    fib = pd.to_numeric(labels["fibrosis_stage"], errors="coerce")
    event = (fib >= 3).astype(float)
    event[fib.isna() | (fib < 0)] = np.nan
    cps = labels["cps"].astype(float)

    # Valid: need event, pseudotime, and CPS
    valid = event.notna() & pseudotime.notna() & cps.notna()
    n_valid = valid.sum()
    log.info("Samples with valid event + pseudotime + CPS: %d", n_valid)

    if n_valid < 40:
        log.warning("Too few samples (%d < 40). Skipping survival curves.", n_valid)
        return pd.DataFrame()

    ev = event[valid].values
    pt = pseudotime[valid].values
    cp = cps[valid].values

    # CPS quartiles
    quartile_edges = np.nanpercentile(cp, [0, 25, 50, 75, 100])
    quartile_labels = CPS_QUARTILE_LABELS
    quartile_idx = np.digitize(cp, quartile_edges[1:-1])  # 0,1,2,3

    # Pseudotime bins
    pt_bins = np.linspace(np.nanmin(pt), np.nanmax(pt), N_PSEUDOTIME_BINS + 1)
    pt_bin_centers = 0.5 * (pt_bins[:-1] + pt_bins[1:])

    curve_rows = []
    for q in range(4):
        q_mask = quartile_idx == q
        n_q = q_mask.sum()
        if n_q < 5:
            continue

        ev_q = ev[q_mask]
        pt_q = pt[q_mask]

        # At each pseudotime bin: fraction of samples at or before this bin
        # that do NOT have F>=3 (i.e., "survival" fraction)
        for b in range(N_PSEUDOTIME_BINS):
            # Samples with pseudotime <= bin upper edge
            at_risk = pt_q <= pt_bins[b + 1]
            n_at_risk = at_risk.sum()
            if n_at_risk == 0:
                continue
            n_event = ev_q[at_risk].sum()
            survival_frac = 1.0 - (n_event / n_at_risk)

            curve_rows.append({
                "cps_quartile": quartile_labels[q],
                "pseudotime_bin": pt_bin_centers[b],
                "pseudotime_upper": pt_bins[b + 1],
                "n_at_risk": n_at_risk,
                "n_event": int(n_event),
                "survival_fraction": survival_frac,
            })

    curves_df = pd.DataFrame(curve_rows)

    if len(curves_df) > 0:
        # Log summary per quartile at final bin
        for q_label in quartile_labels:
            q_data = curves_df[curves_df["cps_quartile"] == q_label]
            if len(q_data) > 0:
                last = q_data.iloc[-1]
                log.info("  %s: final survival=%.3f, n_at_risk=%d, n_events=%d",
                         q_label, last["survival_fraction"],
                         int(last["n_at_risk"]), int(last["n_event"]))

    # Log-rank-like test: compare event rates across CPS quartiles
    # Using chi-squared test on event proportions per quartile
    contingency_table = np.zeros((4, 2))  # quartile x (no_event, event)
    for q in range(4):
        q_mask = quartile_idx == q
        n_q = q_mask.sum()
        if n_q > 0:
            n_event_q = ev[q_mask].sum()
            contingency_table[q, 0] = n_q - n_event_q
            contingency_table[q, 1] = n_event_q

    # Remove empty rows
    nonzero_rows = contingency_table.sum(axis=1) > 0
    if nonzero_rows.sum() >= 2:
        chi2, p_val, dof, expected = stats.chi2_contingency(
            contingency_table[nonzero_rows]
        )
        log.info("Log-rank proxy (chi2 on event proportions): chi2=%.2f, p=%.2e, dof=%d",
                 chi2, p_val, dof)
        # Append this as a metadata row
        curves_df.attrs["logrank_chi2"] = chi2
        curves_df.attrs["logrank_pval"] = p_val

    return curves_df


# ══════════════════════════════════════════════════════════════════════════════
# PART 3: Pseudo-C-index
# ══════════════════════════════════════════════════════════════════════════════
def compute_c_index(features, labels, pseudotime):
    """
    Concordance index: does the model's predicted risk agree with
    the ordering of events along pseudotime?

    Time = pseudotime, Event = F>=3, Predicted risk = logistic regression score.
    """
    log.info("=" * 70)
    log.info("PART 3: Pseudo-C-index (Concordance)")
    log.info("=" * 70)

    fib = pd.to_numeric(labels["fibrosis_stage"], errors="coerce")
    event = (fib >= 3).astype(float)
    event[fib.isna() | (fib < 0)] = np.nan

    # Feature set
    feature_cols = [
        "ct_Stellate", "ct_stellate_gt4pct", "ct_fibrogenic_axis",
        "ct_immune_infiltration", "clin_sex",
        "trans_F2_to_F3", "trans_F3_to_F4",
        "hcc_genetic_score",
    ]
    # Add top 10 div genes
    div_cols = [c for c in features.columns if c.startswith("div_")][:10]
    feature_cols.extend(div_cols)
    feature_cols = [c for c in feature_cols if c in features.columns]

    valid = event.notna() & pseudotime.notna()
    for c in feature_cols:
        valid = valid & features[c].notna()

    n_valid = valid.sum()
    log.info("Samples with valid event + pseudotime + all features: %d", n_valid)

    if n_valid < 50:
        log.warning("Too few samples (%d < 50). Skipping C-index.", n_valid)
        return pd.DataFrame()

    X = features.loc[valid, feature_cols].values.astype(float)
    y = event[valid].values
    pt = pseudotime[valid].values

    # Impute NaN (should be minimal after validity check, but safety)
    for j in range(X.shape[1]):
        nan_mask = np.isnan(X[:, j])
        if nan_mask.any():
            X[nan_mask, j] = np.nanmedian(X[:, j])

    # Fit logistic regression to get predicted risk
    scaler = StandardScaler()
    X_s = scaler.fit_transform(X)
    lr = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=5000,
                            random_state=SEED)
    lr.fit(X_s, y)
    predicted_risk = lr.predict_proba(X_s)[:, 1]

    # Compute C-index
    log.info("Computing concordance index ...")
    c_idx, n_conc, n_disc, n_tied, n_pairs = concordance_index_fast(
        pt, predicted_risk, y.astype(bool)
    )

    log.info("Pseudo-C-index: %.4f", c_idx)
    log.info("  Concordant: %d, Discordant: %d, Tied: %d, Total pairs: %d",
             n_conc, n_disc, n_tied, n_pairs)

    # Also compute C-index for individual features (univariate)
    c_rows = []
    c_rows.append({
        "model": "multivariate_logistic",
        "features": "; ".join(feature_cols),
        "c_index": c_idx,
        "n_concordant": n_conc,
        "n_discordant": n_disc,
        "n_tied": n_tied,
        "n_pairs": n_pairs,
        "n_samples": n_valid,
        "note": "Pseudo-C-index: pseudotime as time, F>=3 as event. NOT true survival.",
    })

    # CPS alone as predictor
    cps = labels.loc[valid, "cps"].astype(float) if "cps" in labels.columns else None
    if cps is not None and cps.notna().sum() > 50:
        cps_valid = cps.notna()
        if cps_valid.sum() > 50:
            c_cps, nc, nd, nt, np_ = concordance_index_fast(
                pt[cps_valid], cps[cps_valid].values, y[cps_valid].astype(bool)
            )
            c_rows.append({
                "model": "CPS_alone",
                "features": "cps",
                "c_index": c_cps,
                "n_concordant": nc,
                "n_discordant": nd,
                "n_tied": nt,
                "n_pairs": np_,
                "n_samples": cps_valid.sum(),
                "note": "CPS composite score alone as risk predictor.",
            })
            log.info("CPS-alone C-index: %.4f (n=%d)", c_cps, cps_valid.sum())

    # Stellate fraction alone
    if "ct_Stellate" in features.columns:
        st = features.loc[valid, "ct_Stellate"]
        st_valid = st.notna()
        if st_valid.sum() > 50:
            c_st, nc, nd, nt, np_ = concordance_index_fast(
                pt[st_valid], st[st_valid].values, y[st_valid].astype(bool)
            )
            c_rows.append({
                "model": "stellate_alone",
                "features": "ct_Stellate",
                "c_index": c_st,
                "n_concordant": nc,
                "n_discordant": nd,
                "n_tied": nt,
                "n_pairs": np_,
                "n_samples": st_valid.sum(),
                "note": "Stellate cell fraction alone as risk predictor.",
            })
            log.info("Stellate-alone C-index: %.4f (n=%d)", c_st, st_valid.sum())

    # HCC score alone
    if "hcc_genetic_score" in features.columns:
        hcc = features.loc[valid, "hcc_genetic_score"]
        hcc_valid = hcc.notna()
        if hcc_valid.sum() > 50:
            c_hcc, nc, nd, nt, np_ = concordance_index_fast(
                pt[hcc_valid], hcc[hcc_valid].values, y[hcc_valid].astype(bool)
            )
            c_rows.append({
                "model": "hcc_score_alone",
                "features": "hcc_genetic_score",
                "c_index": c_hcc,
                "n_concordant": nc,
                "n_discordant": nd,
                "n_tied": nt,
                "n_pairs": np_,
                "n_samples": hcc_valid.sum(),
                "note": "HCC molecular score alone as risk predictor.",
            })
            log.info("HCC-score-alone C-index: %.4f (n=%d)", c_hcc, hcc_valid.sum())

    c_df = pd.DataFrame(c_rows)
    return c_df


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    log.info("=" * 70)
    log.info("Script 153: Pseudo-Survival on Pseudotime")
    log.info("=" * 70)
    log.info("")
    log.info("CRITICAL CAVEAT:")
    log.info("  Pseudotime != calendar time. This is molecular trajectory analysis,")
    log.info("  NOT true survival prediction. Frame results as: 'features associated")
    log.info("  with advanced disease position in molecular trajectory space.'")
    log.info("  True survival analysis requires longitudinal cohorts.")
    log.info("")

    # Load data
    features = safe_load(FEATURE_MATRIX, index_col="sample_id")
    labels = safe_load(PSEUDO_LABELS, index_col="sample_id")

    if features is None or labels is None:
        log.error("Required input files missing. Exiting.")
        sys.exit(1)

    log.info("Features: %s, Labels: %s", features.shape, labels.shape)

    # Load consensus pseudotime
    cpt = safe_load(CONSENSUS_PT)
    if cpt is not None:
        pseudotime = cpt.set_index("sample_id")["pseudotime_consensus"].reindex(features.index)
        log.info("Pseudotime: %d non-NaN / %d total", pseudotime.notna().sum(), len(pseudotime))
    else:
        log.error("Consensus pseudotime not found. Exiting.")
        sys.exit(1)

    # ── Part 1: Hazard ratios ─────────────────────────────────────────────
    hr_df = compute_hazard_ratios(features, labels, pseudotime)
    if len(hr_df) > 0:
        hr_path = os.path.join(OUT_DIR, "survival_hazard_ratios.csv")
        hr_df.to_csv(hr_path, index=False)
        log.info("Saved: %s (%d rows)", os.path.basename(hr_path), len(hr_df))

    # ── Part 2: Survival curves by CPS quartile ──────────────────────────
    curves_df = survival_curves_by_cps(labels, pseudotime)
    if len(curves_df) > 0:
        curves_path = os.path.join(OUT_DIR, "survival_curves_by_cps.csv")
        curves_df.to_csv(curves_path, index=False)
        log.info("Saved: %s (%d rows)", os.path.basename(curves_path), len(curves_df))

    # ── Part 3: Pseudo-C-index ────────────────────────────────────────────
    c_df = compute_c_index(features, labels, pseudotime)
    if len(c_df) > 0:
        c_path = os.path.join(OUT_DIR, "pseudo_survival_c_index.csv")
        c_df.to_csv(c_path, index=False)
        log.info("Saved: %s (%d rows)", os.path.basename(c_path), len(c_df))

    elapsed = time.time() - t0
    log.info("")
    log.info("Done in %.1f seconds.", elapsed)


if __name__ == "__main__":
    main()
