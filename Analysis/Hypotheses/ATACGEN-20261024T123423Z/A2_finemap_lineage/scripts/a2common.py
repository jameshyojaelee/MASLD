"""Shared helpers for A2: bins, peak classes, matched expectation, block bootstrap."""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
A2 = ROOT / "Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage"
W = A2 / "work"
RES = Path(os.environ.get("A2_RES", A2 / "results"))
LIN = ["hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"]
NONHEP = LIN[1:]
SEED = 20261024
NBOOT = 2000
MIN_BIN = 20
TSS_EDGES = [1_000, 5_000, 20_000, 100_000]
MAF_EDGES = [0.01, 0.05, 0.20]


def peak_classes(df, prefix="in_"):
    """Add int8 class indicator columns from boolean in-peak flags (prefix in_ or inb_)."""
    tag = "" if prefix == "in_" else "b"
    m = np.stack([df[f"{prefix}{l}"].to_numpy(bool) for l in LIN], axis=1)
    k = m.sum(axis=1)
    out = {}
    for j, l in enumerate(LIN):
        out[f"any{tag}_{l}"] = m[:, j]
        out[f"spec{tag}_{l}"] = m[:, j] & (k == 1)
        out[f"shared{tag}_{l}"] = m[:, j] & (k >= 2)
    h = m[:, 0]
    n = m[:, 1:].any(axis=1)
    out[f"Honly{tag}"] = h & ~n
    out[f"Nonly{tag}"] = ~h & n
    out[f"HN{tag}"] = h & n
    out[f"anypeak{tag}"] = k >= 1
    for c, v in out.items():
        df[c] = v.astype(np.int8)
    return list(out)


def add_bins(df, ld_edges):
    """TSS x MAF x LD-tertile bins. ld_edges: {ancestry: (q1, q2)} from window variants."""
    tss = np.digitize(df.tss_dist.to_numpy(), TSS_EDGES)
    maf = df.maf.to_numpy(float)
    mb = np.where(np.isnan(maf), 4, np.digitize(np.nan_to_num(maf, nan=0.0), MAF_EDGES))
    ld = df.ld.to_numpy(float)
    lb = np.full(len(df), 3, np.int64)
    for anc, (q1, q2) in ld_edges.items():
        sel = (df.ancestry.to_numpy() == anc) & ~np.isnan(ld)
        lb[sel] = np.digitize(ld[sel], [q1, q2])
    df["tss_bin"] = tss.astype(np.int8)
    df["maf_bin"] = mb.astype(np.int8)
    df["ld_bin"] = lb.astype(np.int8)
    df["bin_full"] = (tss * 20 + mb * 4 + lb).astype(np.int16)
    df["bin_tm"] = (tss * 5 + mb).astype(np.int16)
    df["bin_t"] = tss.astype(np.int16)
    return df


def ld_tertiles(df):
    return {a: tuple(np.nanquantile(s.ld.to_numpy(float), [1 / 3, 2 / 3]))
            for a, s in df.groupby("ancestry") if s.ld.notna().any()}


def matched_expectation(win, weights, key, classes, levels=("bin_full", "bin_tm", "bin_t")):
    """Expected class probability per weighted variant, from same-window same-bin variants.

    win: window variants (one row per variant per window) with `key`, bin columns, classes.
    weights: rows of win to score (subset, same columns) with column `w`.
    Returns weights with exp_<class> columns and `level_used`.
    """
    stats = {}
    for lev in list(levels) + ["window"]:
        gk = [key] if lev == "window" else [key, lev]
        agg = win.groupby(gk, observed=True)[classes].sum()
        agg["n"] = win.groupby(gk, observed=True).size()
        stats[lev] = agg
    out = weights.copy()
    chosen = np.full(len(out), "", object)
    exp = np.full((len(out), len(classes)), np.nan)
    for lev in list(levels) + ["window"]:
        gk = [key] if lev == "window" else [key, lev]
        st = stats[lev]
        j = out[gk].merge(st, left_on=gk, right_index=True, how="left")
        n = j["n"].to_numpy(float)
        take = (chosen == "") & ((n >= MIN_BIN) | (lev == "window"))
        exp[take] = j.loc[take, classes].to_numpy(float) / n[take, None]
        chosen[take] = lev
    for i, c in enumerate(classes):
        out[f"exp_{c}"] = exp[:, i]
    out["level_used"] = chosen
    return out


def unit_summary(df, unit, classes, prefix_o="obs_", prefix_e="exp_"):
    """Mean O and E per inferential unit (block or gene)."""
    cols = [prefix_o + c for c in classes] + [prefix_e + c for c in classes]
    return df.groupby(unit)[cols].mean()


def boot_ratio(units, classes, seed=SEED, nboot=NBOOT, prefix_o="obs_", prefix_e="exp_"):
    """Ratio sum O / sum E over units with unit bootstrap. Returns dict of arrays."""
    O = units[[prefix_o + c for c in classes]].to_numpy(float)
    E = units[[prefix_e + c for c in classes]].to_numpy(float)
    n = len(units)
    rng = np.random.default_rng(seed)
    C = np.zeros((nboot, n))
    for b in range(nboot):
        C[b] = np.bincount(rng.integers(0, n, n), minlength=n)
    so, se = C @ O, C @ E
    with np.errstate(divide="ignore", invalid="ignore"):
        rb = so / se
    return {"obs": O.mean(0), "exp": E.mean(0), "ratio": O.sum(0) / E.sum(0), "boot": rb}


def ci(x, lo=0.025, hi=0.975):
    """Percentile interval of bootstrap ratios. Drops 0/0 draws (NaN) but keeps 0 and inf,
    and uses observed order statistics so an infinite tail is reported as inf."""
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return np.nan, np.nan
    return (float(np.quantile(x, lo, method="inverted_cdf")),
            float(np.quantile(x, hi, method="inverted_cdf")))


def boot_p(logr):
    """Two-sided bootstrap p for log ratio = 0 (2 x smaller tail, +1 corrected).
    Draws with log ratio -inf or +inf count in their tail; only 0/0 draws are dropped."""
    logr = np.asarray(logr, float)
    logr = logr[~np.isnan(logr)]
    b = len(logr)
    if b == 0:
        return np.nan
    lo = (np.sum(logr <= 0) + 1) / (b + 1)
    hi = (np.sum(logr >= 0) + 1) / (b + 1)
    return min(1.0, 2 * min(lo, hi))


def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    r = pv[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(n)
    qq[o] = np.minimum(r, 1)
    q[ok] = qq
    return q
