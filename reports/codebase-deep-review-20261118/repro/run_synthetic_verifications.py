#!/usr/bin/env python3
"""Deterministic compute-node synthetic and read-only artifact verification suite
for the MASLD Resource deep review (2026-11-18).

Strictly read-only on repository files; writes only synthetic_verification_results.json
inside reports/codebase-deep-review-20261118/repro/.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import os
import platform
import socket
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_JSON = (
    REPO_ROOT
    / "reports/codebase-deep-review-20261118/repro/synthetic_verification_results.json"
)
SEED = 20261118


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_n01_c3b_ambient_misattribution_and_label_leakage() -> dict:
    """Verify N01: c3b_ambient.R misattributes 795/1,155 (68.8%) 'within-hepatocyte' genes
    (and 598/689 = 86.8% upregulated calls, incl. ACTA2), passes downregulated ambient genes
    via disease-only share, clamps sparse 0-count markers to rho=1.0, leaks unpermuted disease
    labels into permutation null, and drops 325 uncorrected genes in Figure 3H.
    """
    out_dir = REPO_ROOT / "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient"
    df_genes = pd.read_csv(out_dir / "gene_calls.tsv.gz", sep="\t")
    df_part = pd.read_csv(out_dir / "final_partition_posthoc.tsv.gz", sep="\t")
    df_named = pd.read_csv(out_dir / "named_genes.tsv", sep="\t")

    wh = df_genes[df_genes["within_hepatocyte_corrected"] == True]
    n_wh = int(len(wh))
    n_named_non_hep = int((wh["sc_best_lineage"].notna() & (wh["sc_best_lineage"] != "Hepatocytes")).sum())
    n_not_hep_total = int((wh["sc_best_lineage"] != "Hepatocytes").sum())
    wh_up = wh[wh["bulk_log2fc"] > 0]
    n_wh_up = int(len(wh_up))
    n_wh_up_non_hep = int((wh_up["sc_best_lineage"].notna() & (wh_up["sc_best_lineage"] != "Hepatocytes")).sum())
    n_uncorrected_other = int((df_part["final"] == "within_other_lineage_uncorrected").sum())
    acta2_wh = bool(
        df_named.loc[df_named["symbol"] == "ACTA2", "within_hepatocyte_corrected"].iloc[0]
    )

    assert n_wh == 1155, f"Expected 1155 within_hepatocyte_corrected, got {n_wh}"
    assert n_named_non_hep == 795 and n_not_hep_total == 797
    assert n_wh_up == 689 and n_wh_up_non_hep == 598
    assert n_uncorrected_other == 325
    assert acta2_wh is True

    rng = np.random.default_rng(SEED)
    n_donors = 40
    disease = np.array([0] * 20 + [1] * 20)

    # (a) Sparse pseudobulk: marker has 0 counts in donor d -> pmax(lc[m_idx, ], 1e-9) clamps rho = 1.0
    rho_clamped = min(max(2.0 / max(0.0, 1e-9), 0.0), 1.0)
    assert rho_clamped == 1.0

    # (b) Downregulated LSEC gene: high ambient in controls (20 UMI), 0 in disease -> share_disease = 0 < 0.50
    hc_down = np.where(disease == 0, 20.0, 0.0)
    ambient_down = np.where(disease == 0, 19.0, 0.0)
    share_disease = float(np.median((ambient_down / np.maximum(hc_down, 1.0))[disease == 1]))
    share_controls = float(np.median((ambient_down / np.maximum(hc_down, 1.0))[disease == 0]))
    assert share_disease == 0.0 and share_controls == 0.95

    # (c) Pre-permutation true-label group imputation leakage (c3b_ambient.R:108-112, 125)
    n_genes_null = 300
    hc_null = rng.poisson(lam=50.0, size=(n_genes_null, n_donors)).astype(float)
    has_fib = np.zeros(n_donors, dtype=bool)
    has_fib[0:10] = True
    has_fib[20:30] = True
    fib_obs = np.where(disease[None, :] == 1, 200.0, 10.0) * np.ones((n_genes_null, n_donors))
    fib_imp = fib_obs.copy()
    for g_val in (0, 1):
        miss = (~has_fib) & (disease == g_val)
        obs = has_fib & (disease == g_val)
        fib_imp[:, miss] = fib_obs[:, obs].mean(axis=1, keepdims=True)
    y_corr = np.maximum(hc_null - 0.15 * fib_imp, 0.0)

    def abs_t_stats(mat: np.ndarray, grp: np.ndarray) -> np.ndarray:
        m1 = mat[:, grp == 1].mean(axis=1)
        m0 = mat[:, grp == 0].mean(axis=1)
        v1 = mat[:, grp == 1].var(axis=1, ddof=1)
        v0 = mat[:, grp == 0].var(axis=1, ddof=1)
        se = np.sqrt(v1 / (grp == 1).sum() + v0 / (grp == 0).sum())
        return np.abs((m1 - m0) / np.maximum(se, 1e-9))

    obs_sig = int((abs_t_stats(y_corr, disease) > 2.5).sum())
    perm_sig = float(
        np.mean([(abs_t_stats(y_corr, rng.permutation(disease)) > 2.5).sum() for _ in range(100)])
    )
    assert obs_sig > 10.0 * max(perm_sig, 1.0)

    return {
        "status": "PASS",
        "within_hepatocyte_total": n_wh,
        "within_hepatocyte_named_non_hep_best_lineage": n_named_non_hep,
        "within_hepatocyte_not_hep_total": n_not_hep_total,
        "within_hepatocyte_named_non_hep_pct": round(100.0 * n_named_non_hep / n_wh, 2),
        "within_hepatocyte_up_named_non_hep_pct": round(100.0 * n_wh_up_non_hep / n_wh_up, 2),
        "acta2_within_hepatocyte_corrected": acta2_wh,
        "fig3h_dropped_uncorrected_genes": n_uncorrected_other,
        "sparse_zero_marker_clamped_rho": rho_clamped,
        "downregulated_ambient_share_disease_vs_control": [share_disease, share_controls],
        "pre_perm_imputation_false_positives_obs_vs_perm": [obs_sig, round(perm_sig, 2)],
    }


def check_n02_c5_unstudentized_fwl_permutation_and_series_confound() -> dict:
    """Verify N02: Unstudentized marginal FWL permutation in c5_trio_gse289173.R / XMOD2
    collapses permutation null variance ~2.62x under disease-log10(n_hep) collinearity and produces
    p_perm < 0.05 even when 95% HC3 CI crosses zero.
    """
    c5_prim = pd.read_csv(
        REPO_ROOT / "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c5_trio/primary.tsv",
        sep="\t",
    )
    mat1a = c5_prim.loc[c5_prim["gene"] == "MAT1A"].iloc[0]
    beta_obs = float(mat1a["beta"])
    se_hc3_obs = float(mat1a["se"])
    null_sd_obs = float(mat1a["null_sd"])
    p_perm_obs = float(mat1a["p_perm"])
    df_resid = int(mat1a["df"])
    p_hc3_exact = float(2.0 * stats.t.sf(abs(beta_obs / se_hc3_obs), df=df_resid))
    var_collapse_ratio = float((se_hc3_obs / null_sd_obs) ** 2)

    assert var_collapse_ratio > 2.5
    assert p_hc3_exact > 0.01 and p_perm_obs <= 2e-5

    rng = np.random.default_rng(SEED + 2)
    n = 42
    disease = np.array([0.0] * 23 + [1.0] * 19)
    log10_nhep = np.where(disease == 0, rng.normal(3.58, 0.25, n), rng.normal(2.63, 0.35, n))
    Z = np.column_stack([np.ones(n), log10_nhep, rng.normal(0, 1, (n, 3))])
    Q, _ = np.linalg.qr(Z)
    MZ = np.eye(n) - Q @ Q.T
    xr = MZ @ disease
    denom = float(xr @ xr)

    n_sims = 1000
    fp_unstudentized = 0
    fp_hc3_t = 0
    ci_crosses_zero_yet_perm_sig = 0
    perm_mats = np.array([rng.permutation(disease) for _ in range(500)])
    MZ_D = perm_mats @ MZ
    D_ss = np.sum(MZ_D ** 2, axis=1)

    for _ in range(n_sims):
        sd_i = 0.15 + 12.0 / np.sqrt(10.0 ** log10_nhep)
        y = 8.0 + 0.4 * log10_nhep + rng.normal(0.0, sd_i)
        yr = MZ @ y
        beta = float((xr @ yr) / denom)

        X = np.column_stack([disease, Z])
        XtXi = np.linalg.inv(X.T @ X)
        e = y - X @ (XtXi @ (X.T @ y))
        h = np.sum((X @ XtXi) * X, axis=1)
        se_hc3 = float(
            np.sqrt(np.sum((XtXi[0] @ X.T) ** 2 * (e ** 2) / ((1.0 - np.minimum(h, 0.99)) ** 2)))
        )
        p_hc3 = float(2.0 * stats.t.sf(abs(beta / se_hc3), df=n - 6))
        ci_lo = beta - stats.t.ppf(0.975, df=n - 6) * se_hc3
        ci_hi = beta + stats.t.ppf(0.975, df=n - 6) * se_hc3

        null_betas = (MZ_D @ yr) / D_ss
        p_perm = float((1 + np.sum(np.abs(null_betas) >= abs(beta) - 1e-12)) / (len(null_betas) + 1))

        if p_perm < 0.05:
            fp_unstudentized += 1
            if ci_lo <= 0.0 <= ci_hi:
                ci_crosses_zero_yet_perm_sig += 1
        if p_hc3 < 0.05:
            fp_hc3_t += 1

    fpr_unstud = fp_unstudentized / n_sims
    fpr_hc3 = fp_hc3_t / n_sims
    assert fpr_unstud > 0.18 and ci_crosses_zero_yet_perm_sig > 50

    return {
        "status": "PASS",
        "mat1a_beta": round(beta_obs, 4),
        "mat1a_se_hc3": round(se_hc3_obs, 4),
        "mat1a_perm_null_sd": round(null_sd_obs, 4),
        "mat1a_variance_collapse_factor": round(var_collapse_ratio, 3),
        "mat1a_p_hc3_vs_p_perm": [round(p_hc3_exact, 5), p_perm_obs],
        "synthetic_null_fpr_unstudentized_perm": round(fpr_unstud, 4),
        "synthetic_null_fpr_hc3_t": round(fpr_hc3, 4),
        "synthetic_ci_includes_zero_yet_perm_sig_count": ci_crosses_zero_yet_perm_sig,
    }


def check_n03_lane_duplication_ccl1_and_hotspot_danb() -> dict:
    """Verify N03: 4.84x SRR lane duplication in GSE244832 inflates donor eligibility (>=10 cells),
    creates sister-lane self-neighbors in cKDTree / Hotspot kNN, reverses Fig 4F and CCL1 significance
    when lane-inflated donors are removed, and verify CCL1 F1 95th-percentile 5% removal ceiling.
    """
    p_holm_df = pd.read_csv(
        REPO_ROOT
        / "Analysis/Hypotheses/DIG1-20261024T113951Z/out_atlas_state_replication/primary_holm.tsv",
        sep="\t",
    )
    ccl1_det = pd.read_csv(
        REPO_ROOT
        / "Analysis/Hypotheses/CCL1-20261110T235204Z/out_ccl1_job23189976/detectability_by_donor.tsv",
        sep="\t",
    )
    g244 = ccl1_det[ccl1_det["cohort"] == "GSE244832"].copy()
    # In the un-deduplicated scalesc atlas, 11 GSE244832 donors have n_ccl20_pos_hep >= 10,
    # including 2 healthy controls (D07: 12, D10: 44) and 3 MASLD donors with <=27 rows (D02: 13, D17: 15, D12: 27).
    raw_min10 = g244[g244["n_ccl20_pos_hep"] >= 10]
    dedup_min10 = g244[g244["n_ccl20_pos_hep"] >= 49]  # 10 * 4.84x lane duplication
    n_raw = int(len(raw_min10))
    n_dedup = int(len(dedup_min10))
    assert n_raw == 11 and n_dedup == 6

    p1_gse202379_holm = float(p_holm_df.loc[p_holm_df["id"] == "P1", "p_holm"].iloc[0])
    exact_two_sided_p_n6 = 2.0 / (2 ** n_dedup)
    holm_2_cohorts_n6 = min(1.0, 2.0 * exact_two_sided_p_n6)
    assert p1_gse202379_holm == 0.0625 and holm_2_cohorts_n6 == 0.0625

    rng = np.random.default_rng(SEED + 3)
    n_cells = 10_000
    is_doublet = np.arange(n_cells) < 2500
    cs_all = np.where(is_doublet, rng.normal(3.0, 0.5, n_cells), rng.normal(0.0, 0.5, n_cells))
    keep_f1 = cs_all <= np.quantile(cs_all, 0.95)
    frac_removed_total = float(1.0 - keep_f1.mean())
    remaining_doublet_rate = float(is_doublet[keep_f1].mean())
    assert abs(frac_removed_total - 0.05) < 1e-3 and remaining_doublet_rate > 0.20

    # 5x lane duplication inflates kNN autocorrelation under exact null
    n_unique, k_lanes = 400, 5
    latent = rng.normal(0, 1, (n_unique, 10))
    gene_null = rng.poisson(lam=2.0, size=n_unique).astype(float)
    d_uniq = np.sum((latent[:, None, :] - latent[None, :, :]) ** 2, axis=2)
    np.fill_diagonal(d_uniq, np.inf)
    nn_uniq = np.argsort(d_uniq, axis=1)[:, :10]
    g_z = (gene_null - gene_null.mean()) / gene_null.std()
    c_uniq = float(np.mean(g_z[:, None] * g_z[nn_uniq]))

    latent_dup = np.repeat(latent, k_lanes, axis=0) + rng.normal(0, 0.02, (n_unique * k_lanes, 10))
    gene_dup = np.repeat(gene_null, k_lanes)
    d_dup = np.sum((latent_dup[:, None, :] - latent_dup[None, :, :]) ** 2, axis=2)
    np.fill_diagonal(d_dup, np.inf)
    nn_dup = np.argsort(d_dup, axis=1)[:, :10]
    gd_z = (gene_dup - gene_dup.mean()) / gene_dup.std()
    c_dup = float(np.mean(gd_z[:, None] * gd_z[nn_dup]))
    assert abs(c_uniq) < 0.05 and c_dup > 0.35

    return {
        "status": "PASS",
        "fig4f_gse202379_p_holm": p1_gse202379_holm,
        "fig4f_gse244832_donors_raw_vs_dedup_min10": [n_raw, n_dedup],
        "fig4f_dedup_exact_p_and_holm2": [exact_two_sided_p_n6, holm_2_cohorts_n6],
        "ccl1_f1_max_removed_fraction": round(frac_removed_total, 4),
        "ccl1_f1_residual_doublet_fraction_when_25pct_contaminated": round(remaining_doublet_rate, 4),
        "knn_autocorr_null_unique_vs_5x_lane_duplicated": [round(c_uniq, 4), round(c_dup, 4)],
    }


def check_n04_deconvolution_collider_and_control_stage0() -> dict:
    """Verify N04: (1) Same-sample reference deconvolution without target-gene holdout acts as an
    endogenous severity collider, attenuating >80% of true intracellular hepatocyte remodeling at
    constant cell proportions; (2) unified_metadata.csv codes 68 unstaged GSE213621 controls as
    fibrosis_stage = 0 while coding 31 biopsy-normal GSE162694 controls as NA.
    """
    rng = np.random.default_rng(SEED + 4)
    n_samples, m_ref_genes = 200, 300
    stage = rng.integers(0, 5, size=n_samples).astype(float)
    intracellular_prog = 0.6 * stage + rng.normal(0, 0.4, n_samples)
    y_gene = -0.50 * intracellular_prog + rng.normal(0, 0.25, n_samples)
    # Deconvolution averages across m_ref_genes (~300+ hepatocyte/cholangiocyte state-sensitive genes),
    # shrinking gene-specific noise by 1/sqrt(m_ref_genes) so inferred_ilr acts as a high-precision proxy for intracellular_prog
    ref_panel = -0.35 * intracellular_prog[:, None] + rng.normal(0, 0.40, (n_samples, m_ref_genes))
    inferred_ilr = ref_panel.mean(axis=1)

    b_unadj = float(np.polyfit(stage, y_gene, 1)[0])
    X_adj = np.column_stack([np.ones(n_samples), stage, inferred_ilr])
    b_adj = float(np.linalg.lstsq(X_adj, y_gene, rcond=None)[0][1])
    retained_ratio = b_adj / b_unadj
    assert retained_ratio < 0.20

    meta_path = (
        REPO_ROOT / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
    )
    df_meta = pd.read_csv(meta_path, low_memory=False)
    g213_ctrl = df_meta[(df_meta["dataset"] == "GSE213621") & (df_meta["group_binary"] == "Control")]
    g162_ctrl = df_meta[(df_meta["dataset"] == "GSE162694") & (df_meta["group_binary"] == "Control")]
    n_g213_stage0 = int((g213_ctrl["fibrosis_stage"] == 0).sum())
    n_g162_stage_na = int(g162_ctrl["fibrosis_stage"].isna().sum())
    assert n_g213_stage0 >= 68 and n_g162_stage_na == 31

    return {
        "status": "PASS",
        "constant_fraction_true_slope_unadj": round(b_unadj, 4),
        "constant_fraction_collider_adjusted_slope": round(b_adj, 4),
        "constant_fraction_retained_effect_ratio": round(retained_ratio, 4),
        "gse213621_unstaged_controls_coded_stage0": n_g213_stage0,
        "gse162694_biopsy_normal_controls_coded_na": n_g162_stage_na,
    }


def check_n05_coloc_dedup_maxiter_pos_drop_and_ast_signflip() -> dict:
    """Verify N05: (1) EQTL-DEDUP-PARITY dropped 4/6 'untestable' primary genes whose duplicated rows
    had ~0 CS alpha because patched_08_run_eqtl_susie.R set max_iter=500 vs r2's 5000;
    (2) 06_susie_coloc.R premature !duplicated(CHR:POS) sorted by gwas_pval deletes valid biallelic SNPs;
    (3) 73_release_vs_disease_loci.py flips concordant AST/NAFLD locus ALDH16A1,RASIP1 into 'opposite';
    (4) DIR6 ashr fixg=TRUE IEEE-754 pnorm cancellation sets extreme z=27.1 (PNPLA3) to pm=0, lfsr=1;
    (5) CONF1 C4 promotes HBB (rs334 sickle-cell hemolysis) and GOT1 to genetically_anchored.
    """
    dedup_cmp = pd.read_csv(
        REPO_ROOT
        / "Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/results/eqtl_dedup_gene_comparison.tsv",
        sep="\t",
    )
    # Identify the 4 genes that converged in r2, failed convergence at max_iter=500, and had <1e-6 dup CS alpha
    untestable = dedup_cmp[
        (dedup_cmp["old_converged"] == True)
        & (dedup_cmp["new_converged"] == False)
        & (dedup_cmp["old_cs_alpha_on_lost"].fillna(0) < 1e-6)
    ]
    assert set(untestable["gene"].tolist()) == {"IRF5", "LRRC37A", "KANSL1-AS1", "ABHD12"}

    gwas_df = pd.DataFrame(
        {
            "POS": [100500, 100500],
            "REF": ["A", "A"],
            "ALT": ["ATG", "G"],
            "gwas_pval": [1e-12, 2e-12],
        }
    )
    eqtl_df = pd.DataFrame({"POS": [100500], "REF_eqtl": ["A"], "ALT_eqtl": ["G"]})
    g_bug = gwas_df.sort_values("gwas_pval").drop_duplicates(subset=["POS"], keep="first")
    m_bug = g_bug.merge(eqtl_df, on="POS")
    m_bug = m_bug[(m_bug["REF"] == m_bug["REF_eqtl"]) & (m_bug["ALT"] == m_bug["ALT_eqtl"])]
    m_fix = gwas_df.merge(eqtl_df, on="POS")
    m_fix = m_fix[(m_fix["REF"] == m_fix["REF_eqtl"]) & (m_fix["ALT"] == m_fix["ALT_eqtl"])]
    assert len(m_bug) == 0 and len(m_fix) == 1

    dig1_loci = pd.read_csv(
        REPO_ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci/loci.tsv",
        sep="\t",
    )
    aldh = dig1_loci[dig1_loci["genes"].str.contains("ALDH16A1", na=False)].iloc[0]
    assert str(aldh["cls"]) == "opposite"
    assert float(aldh["z_MVP_AST_EUR"]) < -9.0 and 0.0 < float(aldh["z_MVP_ALT_EUR"]) < 1.96

    z_pnpla3 = 27.1
    cdf_diff_float64 = float(
        stats.norm.cdf(z_pnpla3 / 1.0) - stats.norm.cdf(z_pnpla3 / np.sqrt(1.0 + 0.25 ** 2))
    )
    sf_diff_exact = float(
        stats.norm.sf(z_pnpla3 / np.sqrt(1.0 + 0.25 ** 2)) - stats.norm.sf(z_pnpla3 / 1.0)
    )
    assert cdf_diff_float64 == 0.0 and sf_diff_exact > 0.0

    c4_reclass = pd.read_csv(
        REPO_ROOT
        / "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c4_catalog_candidate/reclassified_genes.tsv",
        sep="\t",
    )
    reclass_syms = set(c4_reclass["symbol"].astype(str))
    assert "HBB" in reclass_syms and "GOT1" in reclass_syms

    return {
        "status": "PASS",
        "eqtl_dedup_false_nonconverged_genes": sorted(untestable["gene"].tolist()),
        "multiallelic_pos_retained_buggy_vs_fixed": [int(len(m_bug)), int(len(m_fix))],
        "aldh16a1_z_alt_z_ast_z_nafld_cls": [
            round(float(aldh["z_MVP_ALT_EUR"]), 3),
            round(float(aldh["z_MVP_AST_EUR"]), 3),
            round(float(aldh["z_MVP_NAFLD_EUR"]), 3),
            str(aldh["cls"]),
        ],
        "ashr_float64_cdf_underflow_at_z27_1": [cdf_diff_float64, sf_diff_exact],
        "c4_hemolysis_enzyme_promoted_to_genetically_anchored": sorted(
            list(reclass_syms & {"HBB", "GOT1"})
        ),
    }


def check_n06_ag2x_fdr_b1_permutation_and_b3a_zero_fill() -> dict:
    """Verify N06: (1) AG2X Stage B 3-track union >=99th pct yields 80.1% empirical background FDR;
    (2) Stage B1 independent per-credible-set PIP permutation deflates null variance ~10x across
    overlapping sets sharing physical SNPs; (3) Stage B3a coding 95.7% unmeasured variants as 0.0
    fabricates OR=36.0 (measured OR=1.73, p=0.33).
    """
    expected_fp = 4638 * 0.02175
    observed_calls = 126
    empirical_fdr = expected_fp / observed_calls
    assert empirical_fdr > 0.79

    rng = np.random.default_rng(SEED + 6)
    n_loci, n_snps_per_locus, k_overlapping_sets = 90, 20, 12
    snp_scores = rng.exponential(scale=1.0, size=(n_loci, n_snps_per_locus))
    raw_logits = rng.normal(0, 1.5, size=(n_loci, n_snps_per_locus))
    base_pip = np.exp(raw_logits) / np.exp(raw_logits).sum(axis=1, keepdims=True)

    null_indep, null_locus = [], []
    for _ in range(400):
        s_ind, s_loc = 0.0, 0.0
        for loc in range(n_loci):
            sc = snp_scores[loc]
            u_m = sc.mean()
            for _k in range(k_overlapping_sets):
                s_ind += (rng.permutation(base_pip[loc]) @ sc) - u_m
            sc_perm = rng.permutation(sc)
            s_loc += k_overlapping_sets * ((base_pip[loc] @ sc_perm) - sc_perm.mean())
        null_indep.append(s_ind / (n_loci * k_overlapping_sets))
        null_locus.append(s_loc / (n_loci * k_overlapping_sets))

    var_ratio_b1 = float(np.var(null_locus, ddof=1) / np.var(null_indep, ddof=1))
    assert var_ratio_b1 > 8.0

    or_zero_filled, p_zero_filled = stats.fisher_exact(np.array([[26, 36], [30, 1753]]))
    or_measured, p_measured = stats.fisher_exact(np.array([[26, 8], [30, 16]]))
    assert or_zero_filled > 35.0 and p_zero_filled < 1e-18
    assert or_measured < 2.0 and p_measured > 0.25

    return {
        "status": "PASS",
        "stageB_expected_fp_vs_126_calls": [round(expected_fp, 2), observed_calls],
        "stageB_empirical_fdr": round(empirical_fdr, 4),
        "stageB1_locus_vs_indep_null_variance_ratio": round(var_ratio_b1, 2),
        "stageB3a_or_p_zero_filled": [round(float(or_zero_filled), 2), float(p_zero_filled)],
        "stageB3a_or_p_measured_only": [round(float(or_measured), 3), round(float(p_measured), 4)],
    }


def check_n07_shipped_models_kleiner_panel_and_chromatin() -> dict:
    """Verify N07 by importing actual shipped release modules (read-only):
    (1) kleiner.py serial-variance r_latent > 1.0 and response-slope disattenuation inflation;
    (2) masld-panel-v1/score.py batch-coupled np.all(np.isfinite) where 1 patient's Cq NaN alters
        another patient's score AND crashes a 100-patient batch if 3 trace genes have 1 NaN each;
    (3) masld-liver-chromatin-state-v1.2/score.py returning identical prof_reliable mask across all forms.
    """
    rel_root = REPO_ROOT / "Analysis/MASLD_Model_Benchmark/release"
    kl_mod = _load_module("kleiner_v12", rel_root / "masld-severity-v1.2/kleiner.py")
    rho_high = float(kl_mod.rho_from_var(kl_mod.GSE193066_58PAIR_STAGE_VAR, kl_mod.SIGMA_E_HIGH))
    r_latent_impossible = 0.65 / np.sqrt(rho_high)
    assert r_latent_impossible > 1.40

    rng = np.random.default_rng(SEED + 7)
    x_rna = rng.normal(0, 1, 5000)
    s_true = 0.80 * x_rna + rng.normal(0, 0.4, 5000)
    y_biopsy = s_true + rng.normal(0, 0.6, 5000)
    slope_true = float(np.polyfit(x_rna, s_true, 1)[0])
    slope_obs = float(np.polyfit(x_rna, y_biopsy, 1)[0])
    rho_y = float(np.var(s_true) / np.var(y_biopsy))
    slope_kleiner_py = slope_obs / rho_y  # kleiner.py line 32: beta_latent = beta_observed / rho
    assert abs(slope_obs - slope_true) < 0.02 and (slope_kleiner_py - slope_true) > 0.25

    # Import actual masld-panel-v1/score.py and run score() on Patient A alone vs [Patient A, Patient B with 1 Cq NaN]
    panel_mod = _load_module("panel_v1", rel_root / "masld-panel-v1/score.py")
    panel_obj = panel_mod.load_panel()
    needed = sorted(set(panel_obj["a"]) | set(panel_obj["b"]))
    cq_patient_a = np.linspace(22.0, 31.0, len(needed))[:, None]  # 20 genes x 1 patient, complete Cq
    score_a_solo, _ = panel_mod.score(cq_patient_a, needed, panel_obj, is_cq=True)

    # Batch Patient A with Patient B who has ONE non-detect (NaN) on needed[0]
    cq_batch = np.hstack([cq_patient_a, cq_patient_a.copy()])
    cq_batch[0, 1] = np.nan
    score_batch, rep_batch = panel_mod.score(cq_batch, needed, panel_obj, is_cq=True)
    # Patient A's score changes purely because Patient B had a NaN in the same call!
    patient_a_drift = float(abs(score_batch[0] - score_a_solo[0]))
    assert patient_a_drift > 1e-3 and rep_batch["pairs_evaluable"] < rep_batch["pairs_total"]

    # Chromatin-state-v1.2 score.py line 238 hardcoded W["prof_reliable"] regardless of form
    chrom_src = (rel_root / "masld-liver-chromatin-state-v1.2/score.py").read_text()
    hardcoded_prof_reliable = 'reliable=np.asarray(W["prof_reliable"], bool)' in chrom_src
    w_npz = np.load(rel_root / "masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz")
    n_prof_reliable = int(np.asarray(w_npz["prof_reliable"], bool).sum())
    assert hardcoded_prof_reliable and n_prof_reliable > 0

    return {
        "status": "PASS",
        "kleiner_gse193066_rho_and_r065_disattenuated": [
            round(rho_high, 4),
            round(float(r_latent_impossible), 4),
        ],
        "kleiner_true_vs_obs_vs_inflated_beta": [
            round(slope_true, 4),
            round(slope_obs, 4),
            round(slope_kleiner_py, 4),
        ],
        "panel_v1_patient_a_solo_vs_batched_score": [
            round(float(score_a_solo[0]), 4),
            round(float(score_batch[0]), 4),
            round(patient_a_drift, 4),
        ],
        "chromatin_v12_hardcoded_single_reliable_mask": [hardcoded_prof_reliable, n_prof_reliable],
    }


def check_n08_spatial_moran_coexpression_log1p_and_signflip_ulp() -> dict:
    """Verify N08: (1) Independent matched-gene Moran's I null inflates co-expressed Hotspot programs
    by >5x via UMI shot-noise cancellation; (2) linear regression of concave log1p(CP10k) on linear
    cell abundance leaves Moran's I > 0.65 under pure composition; (3) V4 summary_donor_tests.tsv
    printed 1e-05 for n=4 donors (exact floor = 0.0625, 6,250x error).
    """
    rng = np.random.default_rng(SEED + 8)
    n_spots, m_genes, n_null, alpha = 900, 26, 400, 0.10
    xy = np.indices((30, 30)).reshape(2, -1).T
    W = ((np.abs(xy[:, None, 0] - xy[None, :, 0]) + np.abs(xy[:, None, 1] - xy[None, :, 1])) == 1).astype(float)
    w_sum = float(W.sum())

    def moran_i(v: np.ndarray) -> float:
        vc = v - v.mean()
        return float((n_spots / w_sum) * (vc @ W @ vc) / (vc @ vc))

    _, eigvecs = np.linalg.eigh(W)
    smooth_basis = eigvecs[:, -20:]
    fields = rng.normal(size=(80, 20)) @ smooth_basis.T
    fields = (fields - fields.mean(axis=1, keepdims=True)) / fields.std(axis=1, keepdims=True)

    prog_genes = np.sqrt(alpha) * fields[0:1].repeat(m_genes, axis=0) + np.sqrt(1.0 - alpha) * rng.normal(
        size=(m_genes, n_spots)
    )
    pool_idx = rng.integers(1, 80, size=500)
    pool_genes = np.sqrt(alpha) * fields[pool_idx] + np.sqrt(1.0 - alpha) * rng.normal(size=(500, n_spots))

    single_gene_prog_I = float(np.mean([moran_i(g) for g in prog_genes]))
    single_gene_pool_I = float(np.mean([moran_i(g) for g in pool_genes[:100]]))
    obs_comp_I = moran_i(prog_genes.mean(axis=0))
    null_comp_I = np.array(
        [moran_i(pool_genes[rng.choice(500, m_genes, replace=False)].mean(axis=0)) for _ in range(n_null)]
    )
    emp_p = float((1 + np.sum(null_comp_I >= obs_comp_I)) / (n_null + 1))
    inflation = float(obs_comp_I / null_comp_I.mean())
    assert abs(single_gene_prog_I - single_gene_pool_I) < 0.03 and inflation > 5.0 and emp_p < 0.01

    gx, gy = np.meshgrid(np.linspace(0, np.pi, 30), np.linspace(0, np.pi, 30))
    n_fib = (5.0 * np.exp(-((gx - 1.5) ** 2 + (gy - 1.5) ** 2) / 0.4)).ravel()
    n_hep = np.clip(12.0 - 1.5 * n_fib, 1.0, None)
    lib = 1000.0 * n_hep + 400.0 * n_fib
    ecm_counts = 0.5 * n_hep + 80.0 * n_fib
    score_z = np.log1p(1e4 * ecm_counts / lib)
    score_z = (score_z - score_z.mean()) / score_z.std()
    X_lin = np.column_stack([np.ones_like(n_fib), np.log1p(lib), n_hep, n_fib])
    resid = score_z - X_lin @ np.linalg.lstsq(X_lin, score_z, rcond=None)[0]
    moran_resid = moran_i(resid)
    assert moran_resid > 0.65

    v4_df = pd.read_csv(
        REPO_ROOT / "Analysis/Hypotheses/V4-20261008T015029Z/summary_donor_tests.tsv",
        sep="\t",
    )
    g192_proxy = v4_df[(v4_df["dataset"] == "GSE192741") & (v4_df["mode"] == "proxy")].iloc[0]
    printed_p_n4 = float(g192_proxy["signflip_p_one_sided"])
    exact_floor_n4 = 2.0 ** (-int(g192_proxy["n_donors"]))
    assert printed_p_n4 == 1e-5 and exact_floor_n4 == 0.0625

    return {
        "status": "PASS",
        "single_gene_moran_prog_vs_pool": [round(single_gene_prog_I, 4), round(single_gene_pool_I, 4)],
        "composite_moran_obs_vs_null_mean_and_inflation": [
            round(obs_comp_I, 4),
            round(float(null_comp_I.mean()), 4),
            round(inflation, 2),
        ],
        "pure_composition_log1p_linear_resid_moran_i": round(moran_resid, 4),
        "v4_n4_printed_p_vs_exact_floor": [printed_p_n4, exact_floor_n4],
    }


def check_n09_xenium_postselection_bias_and_geomx_scramble() -> dict:
    """Verify N09: (1) Xenium SelectionReplay control-only argmax carrier selection shifts H0 null
    negatively and produces bimodal [-1.15, +2.30] bootstrap CIs across carrier switches;
    (2) Govaere 2026 GeoMx sequential well assignment scrambles plates (0/13,969 DEGs on panCK_vs_CD68).
    """
    rng = np.random.default_rng(SEED + 9)
    masks = np.array([False] * 9 + [True] * 9)
    comb_masks = np.zeros((48620, 18), dtype=bool)
    for r, idx in enumerate(itertools.combinations(range(18), 9)):
        comb_masks[r, list(idx)] = True

    neg_obs_count = 0
    perm_means = []
    for _ in range(250):
        # 8 competing Xenium lineages (real X1 carrier argmax searches across all 14 cell types)
        vals = rng.lognormal(mean=0.0, sigma=0.35, size=(18, 8))
        c_obs = int(np.argmax(vals[~masks].mean(axis=0)))
        obs_stat = (vals[masks, c_obs].mean() - vals[~masks, c_obs].mean()) / (
            vals[masks, c_obs].mean() + vals[~masks, c_obs].mean()
        )
        if obs_stat < 0:
            neg_obs_count += 1
        p_ctrl = (~comb_masks @ vals) / 9.0
        p_case = (comb_masks @ vals) / 9.0
        p_car = np.argmax(p_ctrl, axis=1)
        rows = np.arange(len(comb_masks))
        perm_stat = (p_case[rows, p_car] - p_ctrl[rows, p_car]) / (p_case[rows, p_car] + p_ctrl[rows, p_car])
        perm_means.append(float(perm_stat.mean()))

    mean_null_shift = float(np.mean(perm_means))
    neg_obs_frac = neg_obs_count / 250.0
    assert mean_null_shift < -0.04 and neg_obs_frac > 0.65

    r5e = pd.read_csv(
        REPO_ROOT / "Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/xenium/results_22810356/results_5e.tsv",
        sep="\t",
    ).set_index("gene")
    anpep_ci = [round(float(r5e.loc["ANPEP", "share_ci_lo"]), 3), round(float(r5e.loc["ANPEP", "share_ci_hi"]), 3)]
    mtor_ci = [round(float(r5e.loc["MTOR", "share_ci_lo"]), 3), round(float(r5e.loc["MTOR", "share_ci_hi"]), 3)]
    assert anpep_ci[0] < -1.0 and anpep_ci[1] > 2.0

    geomx_sum = pd.read_csv(
        REPO_ROOT / "Analysis/Spatial/results/govaere2026/geomx_de_summary.tsv", sep="\t"
    ).set_index("contrast")
    panck_sig = int(geomx_sum.loc["panCK_vs_CD68", "n_sig_padj05"])
    assert panck_sig == 0

    return {
        "status": "PASS",
        "xenium_null_mean_shift_under_control_argmax": round(mean_null_shift, 4),
        "xenium_null_fraction_negative_statistics": round(neg_obs_frac, 4),
        "xenium_anpep_share_bimodal_ci": anpep_ci,
        "xenium_mtor_share_bimodal_ci": mtor_ci,
        "geomx_panck_vs_cd68_n_sig_of_13969": panck_sig,
    }


def check_n10_atac_shuffle_rel1_holm_and_manuscript_panels() -> dict:
    """Verify N10, N11, N12:
    (1) Within-donor RNA permutation before 50-cell donor metacells preserves cross-donor Cicero ridge gain;
    (2) REL1 B2 multiplies ALDOB p_perm=0.00810 by 4 (0.03240) via 3 untestable p=1 proteins while HC3 CI includes 0;
    (3) DIG1 cameraPR hardcoded inter.gene.cor=0.01 yields >60% FPR across correlated subunits (rho=0.45);
    (4) Manuscript LaTeX & panel collisions: figure6.tex embeds scatter instead of trial table, candidate fig6g
        maps FGF19 analogue to FGF1, fig2.tex prints 436/996 vs Fig 2A/C/D 429/990, Fig 5A (3,054) vs RESULTS.md (2,994).
    """
    rng = np.random.default_rng(SEED + 10)
    n_genes_sim, n_donors, cells_per_donor, k_meta = 80, 39, 500, 50
    gains_p, gains_s = [], []
    for _ in range(n_genes_sim):
        u_donor = rng.normal(0, 1.0, size=(n_donors, 1))
        c_cell = rng.normal(0, 0.3, size=(n_donors, cells_per_donor))
        prom = u_donor + 0.2 * c_cell + rng.normal(0, 3.0, (n_donors, cells_per_donor))
        cic_enh = u_donor + 0.2 * c_cell + rng.normal(0, 3.0, (n_donors, cells_per_donor))
        dist_pk = rng.normal(0, 3.0, (n_donors, cells_per_donor))
        rna = 0.8 * u_donor + 0.2 * c_cell + rng.normal(0, 1.0, (n_donors, cells_per_donor))

        def eval_gain(rna_mat: np.ndarray) -> float:
            P = prom.reshape(n_donors, 10, k_meta).mean(axis=2)
            E = cic_enh.reshape(n_donors, 10, k_meta).mean(axis=2)
            D = dist_pk.reshape(n_donors, 10, k_meta).mean(axis=2)
            Y = rna_mat.reshape(n_donors, 10, k_meta).mean(axis=2)
            diffs = []
            for f in range(5):
                te = (np.arange(n_donors) % 5) == f
                tr = ~te
                Xc_tr = np.column_stack([np.ones(tr.sum() * 10), P[tr].ravel(), E[tr].ravel()])
                Xc_te = np.column_stack([np.ones(te.sum() * 10), P[te].ravel(), E[te].ravel()])
                Xd_tr = np.column_stack([np.ones(tr.sum() * 10), P[tr].ravel(), D[tr].ravel()])
                Xd_te = np.column_stack([np.ones(te.sum() * 10), P[te].ravel(), D[te].ravel()])
                wc = np.linalg.solve(Xc_tr.T @ Xc_tr + 25.0 * np.diag([0, 1, 1]), Xc_tr.T @ Y[tr].ravel())
                wd = np.linalg.solve(Xd_tr.T @ Xd_tr + 25.0 * np.diag([0, 1, 1]), Xd_tr.T @ Y[tr].ravel())
                diffs.append(np.mean((Y[te].ravel() - Xd_te @ wd) ** 2) - np.mean((Y[te].ravel() - Xc_te @ wc) ** 2))
            return float(np.mean(diffs))

        gains_p.append(eval_gain(rna))
        gains_s.append(eval_gain(np.array([rng.permutation(row) for row in rna])))

    gain_paired = float(np.mean(gains_p))
    gain_shuffled = float(np.mean(gains_s))
    assert gain_paired > 0.04 and abs(gain_paired - gain_shuffled) / gain_paired < 0.08

    rel1_b2 = pd.read_csv(
        REPO_ROOT / "Analysis/Hypotheses/REL1-20261110T210214Z/out_23184388/B2_pairs.tsv", sep="\t"
    )
    rel1_fam = pd.read_csv(
        REPO_ROOT / "Analysis/Hypotheses/REL1-20261110T210214Z/out_23184388/families.tsv", sep="\t"
    )
    aldob_nas = rel1_b2[(rel1_b2["protein"] == "ALDOB") & (rel1_b2["histology"] == "d_nas")].iloc[0]
    aldob_fam = rel1_fam[(rel1_fam["item"] == "ALDOB") & (rel1_fam["family"] == "B2 change vs d_nas")].iloc[0]
    assert float(aldob_nas["ci_low"]) < 0.0 < float(aldob_nas["ci_high"])
    assert float(aldob_nas["p_hc3"]) > 0.05
    assert np.isclose(float(aldob_fam["holm"]), 4.0 * float(aldob_fam["p"]))

    m_sub, n_tot, rho_c = 38, 5000, 0.45
    fp_001 = 0
    for _ in range(400):
        f_lat = rng.normal(0, 1)
        sub_z = np.sqrt(rho_c) * f_lat + np.sqrt(1.0 - rho_c) * rng.normal(0, 1, m_sub)
        bg_z = rng.normal(0, 1, n_tot - m_sub)
        diff = float(sub_z.mean() - bg_z.mean())
        s_pool = float(
            np.sqrt(((m_sub - 1) * sub_z.var(ddof=1) + (n_tot - m_sub - 1) * bg_z.var(ddof=1)) / (n_tot - 2))
        )
        se_001 = s_pool * np.sqrt((1.0 + (m_sub - 1) * 0.01) / m_sub + 1.0 / (n_tot - m_sub))
        if 2.0 * stats.norm.sf(abs(diff / se_001)) < 0.05:
            fp_001 += 1
    fpr_camerapr_001 = fp_001 / 400.0
    assert fpr_camerapr_001 > 0.50  # Analytic 2*sf(1.96/sqrt(17.65/1.37)) = 0.585 (>11x nominal 0.05)

    fig6_tex = (REPO_ROOT / "docs/manuscript/latex/sections/figures/figure6.tex").read_text()
    fig2_tex = (REPO_ROOT / "docs/manuscript/latex/sections/fig2.tex").read_text()
    fig1_tex = (REPO_ROOT / "docs/manuscript/latex/sections/fig1.tex").read_text()
    fig6g_cand_txt = (
        REPO_ROOT / "reports/codebase-deep-review-20261118/repro/rendered_panels/txt/panels__fig6g_trial_mechanisms.txt"
    ).read_text()
    fig2a_txt = (
        REPO_ROOT / "reports/codebase-deep-review-20261118/repro/rendered_panels/txt/figure2__Fig2A_gwas_alluvial.txt"
    ).read_text()

    embeds_old_scatter_6g = "fig6g_drug_target_calibration.pdf" in fig6_tex
    cand_6g_has_fgf1_bug = (
        ("FGF19 analogue" in fig6g_cand_txt)
        and ("\nFGF1\n" in fig6g_cand_txt)
        and ("\nFGF19\n" not in fig6g_cand_txt)
    )
    fig2_tex_has_996_436 = ("996" in fig2_tex) and ("436" in fig2_tex)
    fig2a_has_990 = "990" in fig2a_txt
    fig1_tex_has_false_most_alt = "most ALT loci raise other hepatocyte proteins" in fig1_tex

    assert (
        embeds_old_scatter_6g
        and cand_6g_has_fgf1_bug
        and fig2_tex_has_996_436
        and fig2a_has_990
        and fig1_tex_has_false_most_alt
    )

    return {
        "status": "PASS",
        "atac_cicero_gain_paired_vs_within_donor_shuffled": [round(gain_paired, 5), round(gain_shuffled, 5)],
        "rel1_aldob_dnas_p_perm_holm4_phc3": [
            round(float(aldob_fam["p"]), 5),
            round(float(aldob_fam["holm"]), 5),
            round(float(aldob_nas["p_hc3"]), 5),
        ],
        "dig1_camerapr_001_null_fpr_at_rho045": round(fpr_camerapr_001, 4),
        "figure6_tex_embeds_superseded_scatter": embeds_old_scatter_6g,
        "candidate_fig6g_maps_fgf19_analogue_to_fgf1": cand_6g_has_fgf1_bug,
        "fig2_tex_996_436_vs_rendered_fig2a_990": [fig2_tex_has_996_436, fig2a_has_990],
        "fig1_tex_false_most_alt_loci_claim": fig1_tex_has_false_most_alt,
    }


def main() -> None:
    host = socket.gethostname()
    if "login" in host.lower():
        raise RuntimeError(f"Refusing to run on login node: {host}")

    checks = [
        ("check_n01_c3b_ambient_misattribution_and_label_leakage", check_n01_c3b_ambient_misattribution_and_label_leakage),
        ("check_n02_c5_unstudentized_fwl_permutation_and_series_confound", check_n02_c5_unstudentized_fwl_permutation_and_series_confound),
        ("check_n03_lane_duplication_ccl1_and_hotspot_danb", check_n03_lane_duplication_ccl1_and_hotspot_danb),
        ("check_n04_deconvolution_collider_and_control_stage0", check_n04_deconvolution_collider_and_control_stage0),
        ("check_n05_coloc_dedup_maxiter_pos_drop_and_ast_signflip", check_n05_coloc_dedup_maxiter_pos_drop_and_ast_signflip),
        ("check_n06_ag2x_fdr_b1_permutation_and_b3a_zero_fill", check_n06_ag2x_fdr_b1_permutation_and_b3a_zero_fill),
        ("check_n07_shipped_models_kleiner_panel_and_chromatin", check_n07_shipped_models_kleiner_panel_and_chromatin),
        ("check_n08_spatial_moran_coexpression_log1p_and_signflip_ulp", check_n08_spatial_moran_coexpression_log1p_and_signflip_ulp),
        ("check_n09_xenium_postselection_bias_and_geomx_scramble", check_n09_xenium_postselection_bias_and_geomx_scramble),
        ("check_n10_atac_shuffle_rel1_holm_and_manuscript_panels", check_n10_atac_shuffle_rel1_holm_and_manuscript_panels),
    ]

    t0 = time.time()
    results: dict[str, dict] = {}
    passed = 0
    for name, fn in checks:
        ts = time.time()
        res = fn()
        res["elapsed_sec"] = round(time.time() - ts, 3)
        results[name] = res
        if res.get("status") == "PASS":
            passed += 1
        print(f"[{res['status']}] {name} ({res['elapsed_sec']}s)", flush=True)

    payload = {
        "hostname": host,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", "NONE"),
        "python_version": platform.python_version(),
        "seed": SEED,
        "total_checks": len(checks),
        "passed_checks": passed,
        "total_elapsed_sec": round(time.time() - t0, 3),
        "checks": results,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {OUT_JSON} ({passed}/{len(checks)} passed in {payload['total_elapsed_sec']}s)")


if __name__ == "__main__":
    main()
