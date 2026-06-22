#!/usr/bin/env python
"""
DIRECTIONAL forward-validation re-test (TEST 1b).

Reuses the per-gene SIGNED effects already computed
(saunders_forward_validation_per_gene.csv) — does NOT reload the h5ad.

Logic:
  atlas_disease_logFC_g = bulk disease logFC for gene g (canonical C2 DEG; sign = disease direction)
  directional_reversal_g = -sign(atlas_disease_logFC_g) * signed_disease_axis_effect_g
    POSITIVE  = CRISPRi knockdown moved hepatocytes OPPOSITE to g's disease direction
                (expected for a true driver: knock down a disease-UP gene -> disease score DOWN)
    NEGATIVE  = knockdown moved cells FURTHER along the disease direction

Tests:
  (a) one-sample sign + Wilcoxon that directional_reversal > 0, across all genes with a
      clear atlas direction, and separately for the disease-UP subset (cleanest prediction)
  (b) Spearman(directional_reversal, convergence_score) + 10k perm null, raw AND partial
      (control essentiality_chronos + log n_cells)
  (c) high- vs low-convergence Mann-Whitney on directional_reversal
  (2) module-24 (NRF2 antioxidant) signed eigengene shift: top movers, NRF2/high-conv panel genes

Env: spatial. Run via env binary.
"""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, wilcoxon, binomtest, rankdata

np.random.seed(42)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUTDIR = f"{ROOT}/Analysis/SingleCell/results_gpu_v2/perturb_forward_validation"
PERGENE = f"{OUTDIR}/saunders_forward_validation_per_gene.csv"
DEG  = f"{ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"
STATS = f"{OUTDIR}/forward_validation_stats.json"

def log(*a): print(*a, flush=True)

# ---------- load existing per-gene signed effects (no h5ad reload) ----------
d = pd.read_csv(PERGENE)
log(f"reusing per-gene effects: {len(d)} genes")
# robust to re-runs: only rename if the original column names are still present
d = d.rename(columns={"disease_axis_effect": "signed_effect",
                      "module24_effect": "module24_signed_effect"})
assert "signed_effect" in d.columns, "expected signed_effect after rename"

# ---------- atlas disease DIRECTION from canonical C2 DEG (per human ortholog) ----------
deg = pd.read_csv(DEG, usecols=["symbol", "logFC", "shrunk_logFC", "padj"])
deg = deg.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
# use shrunk_logFC for direction (ashr-shrunk effect; sign identical to logFC where both exist)
d["atlas_disease_logFC"] = d["human_ortholog"].map(deg["shrunk_logFC"])
d["atlas_disease_logFC_raw"] = d["human_ortholog"].map(deg["logFC"])
d["atlas_disease_padj"] = d["human_ortholog"].map(deg["padj"])
# fall back to raw logFC if shrunk missing
d["atlas_disease_logFC"] = d["atlas_disease_logFC"].fillna(d["atlas_disease_logFC_raw"])

# clear direction = a defined, non-tiny logFC (avoid sign noise around 0)
DIR_FLOOR = 0.10  # |logFC| floor to call a confident disease direction
d["has_clear_dir"] = d["atlas_disease_logFC"].abs() >= DIR_FLOOR
d["atlas_dir"] = np.sign(d["atlas_disease_logFC"])

# directional reversal
d["directional_reversal"] = -d["atlas_dir"] * d["signed_effect"]

clear = d[d["has_clear_dir"] & d["directional_reversal"].notna()].copy()
up_sub = clear[clear["atlas_dir"] > 0].copy()    # disease-UP genes: cleanest prediction
dn_sub = clear[clear["atlas_dir"] < 0].copy()
log(f"genes with clear atlas direction (|logFC|>={DIR_FLOOR}): {len(clear)} "
    f"(disease-UP={len(up_sub)}, disease-DOWN={len(dn_sub)})")

# also a strict subset: genes that are themselves significant Tier-1 DEGs
strict = d[(d["atlas_disease_padj"] < 0.05) & (d["atlas_disease_logFC"].abs() > 0.5) &
           d["directional_reversal"].notna()].copy()
strict_up = strict[strict["atlas_dir"] > 0]
log(f"strict significant-DEG perturbed genes: {len(strict)} (UP={len(strict_up)})")

# ---------- (a) one-sample sign + Wilcoxon: directional_reversal > 0 ----------
log("\n===== (a) ONE-SAMPLE: is knockdown directionally REVERSING disease? =====")
def one_sample(x, label):
    x = np.asarray(x.dropna())
    n = len(x)
    if n < 3:
        log(f"  [{label}] n={n} too small"); return dict(n=n)
    npos = int((x > 0).sum()); nneg = int((x < 0).sum())
    med = float(np.median(x)); mean = float(np.mean(x))
    sgn = binomtest(npos, npos + nneg, 0.5, alternative="greater").pvalue if (npos+nneg) else np.nan
    try:
        w_stat, w_p = wilcoxon(x, alternative="greater")
        w_p = float(w_p)
    except Exception:
        w_p = np.nan
    log(f"  [{label}] n={n}  pos/neg={npos}/{nneg}  median={med:+.4f}  mean={mean:+.4f}"
        f"  sign-test p(>0)={sgn:.4f}  Wilcoxon p(>0)={w_p:.4f}")
    return dict(n=n, n_pos=npos, n_neg=nneg, median=med, mean=mean,
                sign_p_greater=float(sgn), wilcoxon_p_greater=w_p)

a_all    = one_sample(clear["directional_reversal"],  "all clear-direction")
a_up     = one_sample(up_sub["directional_reversal"], "disease-UP subset (cleanest)")
a_dn     = one_sample(dn_sub["directional_reversal"], "disease-DOWN subset")
a_strict = one_sample(strict["directional_reversal"], "strict Tier-1 DEG subset")

# ---------- (b) Spearman(directional_reversal, convergence_score) raw + partial ----------
log("\n===== (b) Spearman(directional_reversal, convergence_score) =====")
def perm_spearman(x, y, n=10000, seed=42):
    rng = np.random.default_rng(seed)
    rho = spearmanr(x, y).statistic
    yv = np.asarray(y)
    null = np.array([spearmanr(x, rng.permutation(yv)).statistic for _ in range(n)])
    p = (np.sum(np.abs(null) >= abs(rho)) + 1) / (n + 1)
    return float(rho), float(p)

b = clear.dropna(subset=["convergence_score", "directional_reversal"]).copy()
rho_raw, p_raw = perm_spearman(b["convergence_score"].values,
                               b["directional_reversal"].values, n=10000)
log(f"  RAW rho = {rho_raw:+.4f}  perm p = {p_raw:.4f}  (n={len(b)})")

def partial_spearman(df, x, y, covars, seed=7):
    import numpy.linalg as la
    dd = df.dropna(subset=[x, y] + covars).copy()
    if len(dd) < 10: return np.nan, np.nan, len(dd)
    R = {c: rankdata(dd[c].values) for c in [x, y] + covars}
    Z = np.column_stack([np.ones(len(dd))] + [R[c] for c in covars])
    def resid(v):
        beta, *_ = la.lstsq(Z, v, rcond=None); return v - Z @ beta
    rx, ry = resid(R[x]), resid(R[y])
    rho = np.corrcoef(rx, ry)[0, 1]
    rng = np.random.default_rng(seed)
    null = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(10000)])
    p = (np.sum(np.abs(null) >= abs(rho)) + 1) / 10001
    return float(rho), float(p), len(dd)

rho_p, p_p, n_p = partial_spearman(clear, "directional_reversal", "convergence_score",
                                   ["essentiality_chronos", "log_n_cells"])
log(f"  PARTIAL rho (ctrl essentiality + log n_cells) = {rho_p:+.4f}  perm p = {p_p:.4f}  (n={n_p})")

# ---------- (c) high vs low convergence Mann-Whitney on directional_reversal ----------
log("\n===== (c) high vs low convergence: Mann-Whitney on directional_reversal =====")
hi = clear.loc[clear["conv_group"] == "high", "directional_reversal"].dropna()
lo = clear.loc[clear["conv_group"] == "low",  "directional_reversal"].dropna()
log(f"  high n={len(hi)} median={hi.median():+.4f} | low n={len(lo)} median={lo.median():+.4f}")
if len(hi) > 2 and len(lo) > 2:
    u_g, p_g = mannwhitneyu(hi, lo, alternative="greater")
    u_t, p_t = mannwhitneyu(hi, lo, alternative="two-sided")
    log(f"  MWU high>low one-sided p={p_g:.4f}  | two-sided p={p_t:.4f}")
else:
    p_g = p_t = np.nan

# ---------- (2) module-24 (NRF2 antioxidant) signed shifts ----------
log("\n===== (2) MODULE-24 (NRF2-SQSTM1 antioxidant) signed eigengene shifts =====")
m = d.dropna(subset=["module24_signed_effect"]).copy()
# drop perturbations whose human ortholog is unknown (NaN -> breaks JSON consumers)
m_named = m[m["human_ortholog"].notna() & (m["human_ortholog"].astype(str) != "nan")].copy()
top_up = m_named.sort_values("module24_signed_effect", ascending=False).head(10)
top_dn = m_named.sort_values("module24_signed_effect", ascending=True).head(10)
log("  TOP perturbations RAISING module-24 (antioxidant program UP):")
for _, r in top_up.iterrows():
    log(f"    {r['human_ortholog']:>10} mod24={r['module24_signed_effect']:+.3f} "
        f"conv={r['convergence_score']:.3f} tier={r['tier']} grp={r['conv_group']} n={int(r['n_cells'])}")
log("  TOP perturbations LOWERING module-24:")
for _, r in top_dn.iterrows():
    log(f"    {r['human_ortholog']:>10} mod24={r['module24_signed_effect']:+.3f} "
        f"conv={r['convergence_score']:.3f} tier={r['tier']} grp={r['conv_group']} n={int(r['n_cells'])}")
nrf2 = m[m["human_ortholog"].isin(["NFE2L2", "KEAP1", "SQSTM1", "TXNRD1", "SOD2",
                                   "AKR1B10", "HKDC1", "NQO1", "GCLC", "GCLM"])]
log(f"  NRF2-axis panel genes present: {nrf2['human_ortholog'].tolist() if len(nrf2) else 'NONE (NFE2L2/KEAP1 not in panel)'}")
# correlation: does module-24 shift track convergence directionally?
mc = m.dropna(subset=["convergence_score"])
rho_m24, p_m24 = perm_spearman(mc["convergence_score"].values, mc["module24_signed_effect"].values, n=10000)
log(f"  Spearman(signed module-24 shift, convergence) = {rho_m24:+.4f} perm p={p_m24:.4f} (n={len(mc)})")

# ---------- write extended CSV ----------
ext_cols = ["target_gene", "human_ortholog", "convergence_score", "convergence_rank",
            "tier", "recommendation", "conv_group",
            "signed_effect", "abs_disease_axis_effect",
            "atlas_disease_logFC", "atlas_disease_padj", "atlas_dir", "has_clear_dir",
            "directional_reversal", "module24_signed_effect", "abs_module24_effect",
            "essentiality_chronos", "coloc_best_susie_pp4", "is_conserved",
            "is_progression_driver", "n_cells", "log_n_cells"]
d_out = d[[c for c in ext_cols if c in d.columns]].sort_values("convergence_score", ascending=False)
d_out.to_csv(PERGENE, index=False)
log(f"\nupdated per-gene CSV (now with directional cols): {PERGENE}")

# ---------- update stats JSON ----------
with open(STATS) as f:
    stats = json.load(f)
stats["directional"] = dict(
    dir_floor_logFC=DIR_FLOOR,
    n_clear_direction=int(len(clear)), n_disease_up=int(len(up_sub)),
    n_disease_down=int(len(dn_sub)), n_strict_tier1=int(len(strict)),
    onesample_all=a_all, onesample_disease_up=a_up,
    onesample_disease_down=a_dn, onesample_strict=a_strict,
    spearman_raw_rho=rho_raw, spearman_raw_perm_p=p_raw, spearman_raw_n=int(len(b)),
    spearman_partial_rho=rho_p, spearman_partial_perm_p=p_p, spearman_partial_n=int(n_p),
    mwu_high_n=int(len(hi)), mwu_low_n=int(len(lo)),
    mwu_high_median=float(hi.median()), mwu_low_median=float(lo.median()),
    mwu_p_greater=float(p_g) if not np.isnan(p_g) else None,
    mwu_p_two_sided=float(p_t) if not np.isnan(p_t) else None,
)
stats["module24"] = dict(
    spearman_signed_vs_convergence_rho=rho_m24, perm_p=p_m24,
    top_raisers=top_up[["human_ortholog", "module24_signed_effect", "convergence_score", "tier"]].to_dict("records"),
    top_lowerers=top_dn[["human_ortholog", "module24_signed_effect", "convergence_score", "tier"]].to_dict("records"),
    nrf2_axis_panel_genes=nrf2["human_ortholog"].tolist(),
)
def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, float) and np.isnan(o):
        return None
    return o
with open(STATS, "w") as f:
    json.dump(_clean(stats), f, indent=2)
log(f"updated stats JSON: {STATS}")
log("\n==== DIRECTIONAL DONE ====")
