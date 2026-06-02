"""Shared donor-level pseudobulk utilities for the scATAC pipeline.

Created 2026-05-30 to fix the pervasive PSEUDOREPLICATION identified in the
2026-05-29 code review (findings F024/F028/F033/F038/F048/F075/F079/F150):
several differential tests treated individual CELLS (~40k) as independent
replicates when the true experimental unit is the DONOR (18 donors, ~6 NORMAL
/ 12 MASL+MASH). Collapsing cells -> donors before testing makes the effective
n the number of donors, not cells.

This mirrors the already-correct pattern in 09_chromvar_propagation.py (per-donor
Wilcoxon) and 14_stage_stratified_da.py (per-donor pseudobulk). Consumers:
03_chromvar_motifs.py, 04b_scenic_grn_from_activity.py, 07b_corrected_da_hepatocytes.py,
13_sex_epigenomic_validation.py, 16_loo_donor_validation.py.

API
---
aggregate_cells_to_donors(X, donor_array, agg="mean"|"sum", min_cells=1)
    -> (donors:list, M:ndarray[n_donor, n_feature], n_cells:ndarray[n_donor])
donor_groupwise_test(M, donors, donor_to_group, group_a, group_b,
                     feature_names=None, min_per_group=2)
    -> pandas.DataFrame with per-feature donor-level Mann-Whitney U + BH-FDR.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

try:
    import scipy.sparse as sp
    _HAVE_SP = True
except Exception:  # pragma: no cover
    _HAVE_SP = False


def aggregate_cells_to_donors(X, donor_array, agg: str = "mean", min_cells: int = 1):
    """Collapse a (n_cells, n_features) matrix to (n_donor, n_features).

    Parameters
    ----------
    X : ndarray or scipy.sparse matrix, shape (n_cells, n_features)
    donor_array : array-like (n_cells,) donor identifiers (any hashable)
    agg : "mean" (per-donor mean; for continuous scores e.g. chromVAR/regulon)
          or "sum" (per-donor pseudobulk counts; for DA peak/tile counts)
    min_cells : drop donors contributing fewer than this many cells

    Returns
    -------
    donors : list of donor ids kept, in sorted order (length n_donor)
    M : ndarray (n_donor, n_features), float64
    n_cells : ndarray (n_donor,) cells aggregated per donor
    """
    if agg not in ("mean", "sum"):
        raise ValueError(f"agg must be 'mean' or 'sum', got {agg!r}")
    donor_array = np.asarray(donor_array)
    if donor_array.shape[0] != X.shape[0]:
        raise ValueError(
            f"donor_array length ({donor_array.shape[0]}) != n_cells ({X.shape[0]})"
        )
    is_sparse = _HAVE_SP and sp.issparse(X)
    donors_all = sorted(pd.unique(donor_array).tolist(), key=lambda d: str(d))
    rows, counts, kept = [], [], []
    for d in donors_all:
        idx = np.where(donor_array == d)[0]
        if idx.size < min_cells:
            continue
        sub = X[idx, :]
        s = np.asarray(sub.sum(axis=0)).ravel().astype(np.float64)
        rows.append(s / idx.size if agg == "mean" else s)
        counts.append(int(idx.size))
        kept.append(d)
    if not rows:
        return [], np.zeros((0, X.shape[1]), dtype=np.float64), np.zeros(0, dtype=int)
    return kept, np.vstack(rows), np.asarray(counts, dtype=int)


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaNs preserved)."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def donor_groupwise_test(
    M: np.ndarray,
    donors: list,
    donor_to_group: dict,
    group_a,
    group_b,
    feature_names=None,
    min_per_group: int = 2,
):
    """Per-feature donor-level Mann-Whitney U (the DONOR is the unit).

    Mirrors 09_chromvar_propagation.py: each donor contributes ONE value per
    feature; we compare the two donor groups (e.g. disease vs healthy). Returns
    a tidy DataFrame with mean_a/mean_b, mean difference, U statistic, raw p,
    BH-FDR padj, and the donor counts actually used.
    """
    from scipy.stats import mannwhitneyu

    M = np.asarray(M, dtype=float)
    n_feat = M.shape[1]
    grp = np.array([donor_to_group.get(d, None) for d in donors], dtype=object)
    a_idx = np.where(grp == group_a)[0]
    b_idx = np.where(grp == group_b)[0]
    n_a, n_b = a_idx.size, b_idx.size
    if feature_names is None:
        feature_names = np.arange(n_feat)

    mean_a = M[a_idx, :].mean(axis=0) if n_a else np.full(n_feat, np.nan)
    mean_b = M[b_idx, :].mean(axis=0) if n_b else np.full(n_feat, np.nan)
    pvals = np.full(n_feat, np.nan)
    stats = np.full(n_feat, np.nan)
    if n_a >= min_per_group and n_b >= min_per_group:
        for j in range(n_feat):
            av, bv = M[a_idx, j], M[b_idx, j]
            # need variation in the pooled values for a defined U test
            if np.all(av == av[0]) and np.all(bv == bv[0]) and av[0] == bv[0]:
                continue
            try:
                u, p = mannwhitneyu(av, bv, alternative="two-sided")
                stats[j], pvals[j] = u, p
            except ValueError:
                continue

    return pd.DataFrame(
        {
            "feature": np.asarray(feature_names),
            f"mean_{group_a}": mean_a,
            f"mean_{group_b}": mean_b,
            "mean_diff": mean_a - mean_b,
            "u_stat": stats,
            "pvalue": pvals,
            "padj": _bh_fdr(pvals),
            "n_donors_a": n_a,
            "n_donors_b": n_b,
        }
    )


# ---------------------------------------------------------------------------
# Self-test: run `python utils_pseudobulk.py`. Verifies (a) aggregation reduces
# cells -> donors with correct means/sums, (b) the donor-level test's effective
# n is the donor count, and (c) a feature with a true donor-level shift is
# detected while a feature that only varies cell-to-cell within donors is not.
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rng = np.random.RandomState(0)
    # Realistic cohort: 18 donors, 13 disease / 5 healthy (the true GSE244832
    # split), ~50 cells each. At n=3 vs 3 the Mann-Whitney two-sided p floors at
    # 0.1, so a smaller cohort cannot distinguish the two features; 13 vs 5 gives
    # min p ~ 2/C(18,5) = 2.3e-4, enough resolution to detect a real donor shift.
    n_dis, n_hlt, cells_per = 13, 5, 50
    ids = [f"D{i:02d}" for i in range(n_dis + n_hlt)]
    donor_ids = np.repeat(ids, cells_per)
    cond = {d: ("disease" if i < n_dis else "healthy") for i, d in enumerate(ids)}
    n_cells = donor_ids.size

    # feature 0: real DONOR-level shift (disease donors offset, per-donor)
    # feature 1: NO donor shift, only within-donor cell noise (the pseudorep trap)
    f0 = np.zeros(n_cells)
    f1 = np.zeros(n_cells)
    for i, d in enumerate(ids):
        m = donor_ids == d
        donor_offset = (2.0 if i < n_dis else 0.0) + rng.randn() * 0.2
        f0[m] = donor_offset + rng.randn(m.sum())   # real donor-level signal
        f1[m] = rng.randn(m.sum())                  # pure cell noise, no donor shift
    X = np.column_stack([f0, f1])

    donors, M, nc = aggregate_cells_to_donors(X, donor_ids, agg="mean")
    assert donors == ids, donors
    assert M.shape == (18, 2), M.shape
    assert np.all(nc == cells_per), nc
    assert abs(M[0, 0] - f0[donor_ids == "D00"].mean()) < 1e-9  # mean sanity

    _, Msum, _ = aggregate_cells_to_donors(X, donor_ids, agg="sum")
    assert abs(Msum[0, 0] - f0[donor_ids == "D00"].sum()) < 1e-6  # sum sanity

    res = donor_groupwise_test(M, donors, cond, "disease", "healthy",
                               feature_names=["real_shift", "cell_noise_only"])
    assert res["n_donors_a"].iloc[0] == 13 and res["n_donors_b"].iloc[0] == 5, \
        "effective n must be donors (13 vs 5), not cells"
    p_real = res.loc[res.feature == "real_shift", "pvalue"].iloc[0]
    p_noise = res.loc[res.feature == "cell_noise_only", "pvalue"].iloc[0]
    assert p_real < 0.05, f"real donor-level shift should be detected, p={p_real}"
    assert p_real < p_noise, (p_real, p_noise)
    # min_per_group guard: leave only 1 healthy donor -> no test
    res2 = donor_groupwise_test(
        M, donors, {**cond, "D13": "x", "D14": "x", "D15": "x", "D16": "x"},
        "disease", "healthy")
    assert np.isnan(res2["pvalue"].iloc[0]), "should refuse with <2 donors/group"

    print("utils_pseudobulk self-test PASSED")
    print(f"  donor-level test effective n = {res['n_donors_a'].iloc[0]} vs "
          f"{res['n_donors_b'].iloc[0]} donors (not {n_cells} cells)")
    print(f"  real_shift p={p_real:.3g}  cell_noise p={p_noise:.3g}")
