"""
spatial_stats.py — donor-aware statistics helpers for the spatial pipeline.

Created 2026-06-01 for the spatial rigor pass (see
docs/reviews/2026-06-01-spatial-transcriptomics-code-review.md). The main
review found pervasive spot/cell-level pseudoreplication: the spatial cohort
is ~5 donors (GSE192741: sample_id JBO014/015/018/019/022) but tests treated
the thousands of *spots* as independent observations. These helpers replace
spot-level tests with donor-aware ones:

  * mixed-effects models (statsmodels MixedLM, donor random intercept) when
    enough donor groups exist to estimate a random-effect variance, and
  * per-donor aggregation + a donor-level test when they do not (n=5 is too
    few groups for a stable random effect — fall back rather than pretend).

Every function reports `method` and `n_donors` so callers can surface the
honest n. With ~5 donors most tests are effectively descriptive; that is the
correct consequence of fixing pseudoreplication, not a bug.

Also provides `ensure_lognorm()` for the raw-counts-in-tests findings
(F016/F038/F094): the deconvolved h5ad `.X` is raw integer counts with no
log/CPM layer, so per-gene/L-R tests were depth-confounded.

Env: spatial (statsmodels 0.14.6, scanpy, squidpy). Pure-python, no R.
"""
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from scipy import stats as _ss

# statsmodels is optional at import time so the module still loads where only
# the aggregation helpers are needed.
try:
    import statsmodels.formula.api as smf
    _HAVE_SM = True
except Exception:  # pragma: no cover
    _HAVE_SM = False

# Minimum number of donor groups before a mixed-effects model is attempted.
# With < this many groups the random-effect variance is unstable, so we fall
# back to per-donor aggregation + a fixed-effects donor-level test.
MIN_DONORS_FOR_MIXED = 6


# --------------------------------------------------------------------------- #
# Normalization
# --------------------------------------------------------------------------- #
def ensure_lognorm(adata, target_sum: float = 1e4, layer_out: str = "lognorm",
                   force: bool = False):
    """Guarantee a log1p-CPM expression layer and return it as a layer name.

    Many spatial scripts read ``adata.X`` directly for per-gene / ligand-
    receptor tests, but ``03b_run_cell2location.py`` sets ``adata.X`` to raw
    integer counts (only a ``'counts'`` layer, no log/CPM). Running Wilcoxon /
    ligrec / Spearman on raw counts confounds the test with per-spot library
    size. This builds a normalized layer once and returns its name; callers
    should pass ``layer=ensure_lognorm(adata)`` to their test.

    Returns the layer name (``layer_out``). Idempotent unless ``force``.
    """
    import scanpy as sc

    if layer_out in adata.layers and not force:
        return layer_out

    # Prefer an explicit raw-counts source; else assume .X is counts.
    if "counts" in adata.layers:
        counts = adata.layers["counts"]
    else:
        counts = adata.X
    tmp = adata.copy()
    tmp.X = counts.copy()
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    sc.pp.log1p(tmp)
    adata.layers[layer_out] = tmp.X
    return layer_out


def _get_expr(adata, gene, layer=None):
    """1-D dense expression vector for one gene from a layer (or .X)."""
    idx = adata.var_names.get_loc(gene)
    mat = adata.layers[layer] if layer is not None else adata.X
    col = mat[:, idx]
    return np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()


# --------------------------------------------------------------------------- #
# Per-donor aggregation
# --------------------------------------------------------------------------- #
def pseudobulk_by_donor(adata, donor_col: str, genes=None, layer=None,
                        agg: str = "mean", obs_cols=None) -> pd.DataFrame:
    """Aggregate spot/cell-level expression to one row per donor.

    Returns a DataFrame indexed by donor with one column per gene (aggregated
    with ``agg``), plus any per-donor ``obs_cols`` (taken as the donor's first
    value — intended for donor-constant columns like condition/disease_stage).
    """
    if genes is None:
        genes = list(adata.var_names)
    donors = adata.obs[donor_col].astype(str).values
    out = {}
    for g in genes:
        v = _get_expr(adata, g, layer=layer)
        s = pd.Series(v).groupby(donors).agg(agg)  # agg is "mean"/"median"/"sum"
        out[g] = s
    df = pd.DataFrame(out)
    df.index.name = donor_col
    if obs_cols:
        meta = (adata.obs[[donor_col] + list(obs_cols)]
                .groupby(donor_col).first())
        df = df.join(meta, how="left")
    return df


# --------------------------------------------------------------------------- #
# Two-group comparison (e.g. Healthy vs Steatotic)
# --------------------------------------------------------------------------- #
def compare_two_groups_by_donor(adata, gene, group_col, donor_col,
                                layer=None, min_donors_for_mixed=MIN_DONORS_FOR_MIXED):
    """Donor-aware two-group test for one gene.

    Chooses the statistically appropriate model from the design:

    * If the grouping factor *varies within* donors (repeated conditions per
      donor) AND there are >= ``min_donors_for_mixed`` donors, fits
      ``expr ~ C(group) + (1|donor)`` (MixedLM) — the mixed model genuinely
      uses within-donor information.
    * Otherwise the factor is *between-donor* (e.g. each donor is entirely
      Healthy or Steatotic, as in this n=5 cohort). A donor random intercept
      cannot separate a between-donor factor (MixedLM returns a nan Wald p),
      and the mixed model adds nothing over aggregation — so we aggregate per
      donor (mean) and run a Mann-Whitney U on the donor-level means. This is
      the honest test for the GSE192741 design.

    Returns a dict with method/pval/effect/n_donors.
    """
    expr = _get_expr(adata, gene, layer=layer)
    df = pd.DataFrame({
        "expr": expr,
        "group": adata.obs[group_col].astype(str).values,
        "donor": adata.obs[donor_col].astype(str).values,
    })
    n_donors = df["donor"].nunique()
    groups = sorted(df["group"].unique())
    # Does the factor vary within at least one donor? If not it is between-donor.
    varies_within = bool((df.groupby("donor")["group"].nunique() > 1).any())
    base = {"gene": gene, "n_donors": int(n_donors), "groups": groups,
            "n_groups": len(groups), "between_donor": (not varies_within)}

    if len(groups) != 2:
        return {**base, "method": "skipped", "pval": np.nan,
                "effect": np.nan, "note": "needs exactly 2 groups"}

    # Within-donor factor: a mixed model is the right tool (the random intercept
    # absorbs donor baseline; the within-donor effect is estimable even at ~5
    # donors given many spots/donor). Require >=3 donors and fall back to
    # aggregation if the fit is degenerate / returns a nan Wald p.
    if _HAVE_SM and varies_within and n_donors >= 3:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                md = smf.mixedlm("expr ~ C(group)", df, groups=df["donor"])
                fit = md.fit(reml=False, method="lbfgs")
            coef = [c for c in fit.params.index if c.startswith("C(group)")]
            if coef and np.isfinite(fit.pvalues[coef[0]]):
                c = coef[0]
                return {**base, "method": "mixedlm",
                        "pval": float(fit.pvalues[c]),
                        "effect": float(fit.params[c])}
        except Exception as e:  # convergence / singular — fall through
            base["mixedlm_error"] = str(e)[:120]

    # Between-donor (or too few donors / mixedlm degenerate): per-donor mean,
    # donor-level Mann-Whitney. Correct test for a donor-level grouping factor.
    agg = df.groupby(["donor", "group"])["expr"].mean().reset_index()
    a = agg.loc[agg["group"] == groups[0], "expr"].values
    b = agg.loc[agg["group"] == groups[1], "expr"].values
    if len(a) < 1 or len(b) < 1:
        return {**base, "method": "donor_agg", "pval": np.nan, "effect": np.nan}
    try:
        u, p = _ss.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        p = np.nan
    return {**base, "method": "donor_agg_mwu",
            "pval": float(p), "effect": float(np.mean(b) - np.mean(a)),
            "n_a": int(len(a)), "n_b": int(len(b))}


def compare_two_groups_many_genes(adata, genes, group_col, donor_col, layer=None,
                                  fdr=True):
    """Vectorized wrapper over compare_two_groups_by_donor for a gene list.

    Returns a DataFrame with per-gene method/pval/effect and BH-FDR padj.
    """
    rows = [compare_two_groups_by_donor(adata, g, group_col, donor_col, layer=layer)
            for g in genes]
    res = pd.DataFrame(rows)
    res["padj"] = np.nan  # always present, even if every pval is nan
    if fdr and "pval" in res and res["pval"].notna().any():
        from statsmodels.stats.multitest import multipletests
        mask = res["pval"].notna()
        res.loc[mask, "padj"] = multipletests(res.loc[mask, "pval"], method="fdr_bh")[1]
    return res


# --------------------------------------------------------------------------- #
# Correlation (e.g. zonation score vs expression) with donor blocking
# --------------------------------------------------------------------------- #
def spearman_by_donor(x, y, donor, min_per_donor=20):
    """Per-donor Spearman rho, then aggregate across donors.

    Computes rho within each donor (>= ``min_per_donor`` obs), then tests
    whether the per-donor rho distribution differs from 0 with a Wilcoxon
    signed-rank test. Reports the mean rho and the donor-level p-value instead
    of a spot-level p-value (which is pseudoreplicated).
    """
    df = pd.DataFrame({"x": np.asarray(x, float), "y": np.asarray(y, float),
                       "donor": np.asarray(donor)})
    rhos = {}
    for d, sub in df.groupby("donor"):
        sub = sub.dropna()
        if len(sub) >= min_per_donor and sub["x"].nunique() > 1 and sub["y"].nunique() > 1:
            rhos[d] = _ss.spearmanr(sub["x"], sub["y"]).correlation
    rho_vals = np.array([v for v in rhos.values() if np.isfinite(v)])
    out = {"n_donors": int(len(rho_vals)),
           "mean_rho": float(np.mean(rho_vals)) if len(rho_vals) else np.nan,
           "per_donor_rho": rhos}
    if len(rho_vals) >= 3:
        try:
            out["pval"] = float(_ss.wilcoxon(rho_vals).pvalue)
            out["method"] = "per_donor_wilcoxon"
        except ValueError:
            out["pval"] = np.nan
            out["method"] = "per_donor_wilcoxon_degenerate"
    else:
        out["pval"] = np.nan
        out["method"] = "too_few_donors_descriptive"
    return out


# --------------------------------------------------------------------------- #
# Spatial autocorrelation (Moran's I) per donor
# --------------------------------------------------------------------------- #
def morans_i_by_donor(adata, genes, donor_col, n_perms=999, n_neighs=6,
                      layer=None, seed=42):
    """Per-donor Moran's I, aggregated across donors.

    Builds the spatial-neighbor graph *within each donor slice* (so spots from
    physically distinct slices are never connected — fixes F229), computes
    Moran's I per donor with squidpy, then aggregates: mean I across donors and
    the fraction of donors with FDR-significant autocorrelation. Replaces the
    pooled-spots single-graph Moran's I (F009/F044/F183).

    Returns a DataFrame indexed by gene.
    """
    import scanpy as sc
    import squidpy as sq

    per_donor = []
    for d, idx in adata.obs.groupby(donor_col).groups.items():
        sub = adata[list(idx)].copy()
        if sub.n_obs < n_neighs + 1:
            continue
        if layer is not None:
            sub.X = sub.layers[layer]
        sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=n_neighs)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sq.gr.spatial_autocorr(sub, mode="moran", genes=list(genes),
                                   n_perms=n_perms, seed=seed)
        m = sub.uns["moranI"].copy()
        m["donor"] = d
        per_donor.append(m)
    if not per_donor:
        return pd.DataFrame()
    allm = pd.concat(per_donor)
    pcol = "pval_norm_fdr_bh" if "pval_norm_fdr_bh" in allm.columns else (
        "pval_norm" if "pval_norm" in allm.columns else None)
    g = allm.groupby(allm.index)
    out = pd.DataFrame({
        "morans_i_mean": g["I"].mean(),
        "n_donors": g["I"].size(),
    })
    if pcol:
        out["frac_donors_sig"] = g.apply(lambda x: float((x[pcol] < 0.05).mean()))
    return out


def per_donor_then_aggregate(adata, donor_col, fn, min_obs=20):
    """Generic: run ``fn(sub_adata)->dict`` per donor, return a tidy DataFrame.

    For donor-blocked versions of squidpy nhood_enrichment / ligrec etc.
    (F036/F208/F210): compute the statistic within each donor, then the caller
    aggregates (e.g. mean z-score + a sign test across donors).
    """
    rows = []
    for d, idx in adata.obs.groupby(donor_col).groups.items():
        sub = adata[list(idx)].copy()
        if sub.n_obs < min_obs:
            continue
        r = fn(sub)
        if r is not None:
            r = dict(r)
            r[donor_col] = d
            rows.append(r)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# COMMOT pathway discovery (F078/F080/F081/F082)
# --------------------------------------------------------------------------- #
# COMMOT (CellChat DB) stores results in a MIX of granularities:
#   * obsp 'commot-CellChat-<PATHWAY>'                  -> true pathway-level
#     matrix (e.g. 'commot-CellChat-TGFb'); the segment after the database
#     prefix contains NO further '-'.
#   * obsp 'commot-CellChat-<LIG>-<REC>'                -> individual LR pair
#     (e.g. 'commot-CellChat-HGF-MET'); has a 2nd '-'.
#   * obsp 'commot-CellChat-total-total'               -> global aggregate.
#   * obsm 'commot-CellChat-sum-sender/-receiver'      -> global sum aggregate.
# The old `key.split('-')[2]` parsing treated every ligand gene as a fake
# pathway and ingested 'total'/'sum' as pathways. These helpers restrict to the
# genuine pathway-level matrices and match config names case-insensitively.

_COMMOT_AGGREGATE_TOKENS = {"total", "sum"}


def get_commot_pathways(adata, db: str = "CellChat"):
    """Return the canonical pathway-level COMMOT names for a database.

    Reads only obsp keys of the form ``commot-<db>-<PATHWAY>`` where the segment
    after ``commot-<db>-`` contains no further ``-`` (so LR-pair matrices like
    ``commot-CellChat-HGF-MET`` are excluded) and is not an aggregate token
    (``total``/``sum``). Cross-checks against the canonical pathway list in
    ``adata.uns['commot-<db>-info']['df_ligrec']['pathway']`` when present.

    Returns a sorted list of real pathway names (e.g. ['CCL','HGF','TGFb',...]).
    """
    prefix = f"commot-{db}-"
    pathways = set()
    for key in getattr(adata, "obsp", {}).keys():
        if not key.startswith(prefix):
            continue
        rest = key[len(prefix):]
        if "-" in rest:           # LR pair (lig-rec) -> not a pathway
            continue
        if rest.lower() in _COMMOT_AGGREGATE_TOKENS:
            continue
        pathways.add(rest)

    # Cross-check / restrict to canonical pathways from the COMMOT info table.
    info_key = f"commot-{db}-info"
    canonical = None
    if info_key in getattr(adata, "uns", {}):
        info = adata.uns[info_key]
        df = info.get("df_ligrec") if isinstance(info, dict) else None
        if df is not None and "pathway" in getattr(df, "columns", []):
            canonical = set(df["pathway"].dropna().astype(str))
    if canonical:
        pathways &= canonical
    return sorted(pathways)


def match_key_pathways(key_pathways, available):
    """Case-insensitive (with WNT-family) match of config pathways to available.

    Config ``key_pathways`` (e.g. ['TGFB','PDGF','CCL','WNT','VEGF','HGF']) use
    upper-case generic names, but CellChat names are case-mixed ('TGFb') and the
    WNT pathway is split into family members ('ncWNT', and any 'WNT*'). Returns
    the subset of ``available`` (canonical-cased) that match any requested key.
    """
    avail_lower = {a.lower(): a for a in available}
    matched = []
    for p in key_pathways:
        pl = str(p).lower()
        if pl in avail_lower:
            matched.append(avail_lower[pl])
        elif pl == "wnt":
            # WNT family: ncWNT / WNT5A / WNT11 etc.
            matched.extend(a for a in available if "wnt" in a.lower())
    # de-dup, preserve canonical order
    seen, out = set(), []
    for a in matched:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out
