#!/usr/bin/env python
"""
Forward-validation TEST 1 for the MASLD convergence atlas.

CLAIM: does the atlas convergence rank PREDICT the functional consequence of
perturbing a gene in hepatocytes, in an INDEPENDENT in-vivo CRISPRi screen
(Saunders 2025) that was NOT designed around our atlas?

Method:
 1. Per hepatocyte: score MASLD disease-UP and disease-DOWN human signatures and
    the hep module-24 (NRF2-SQSTM1 antioxidant) program on the mouse-ortholog-mapped
    gene sets (scanpy score_genes = control-binned mean-z; AUCell-style robust scorer).
    disease-axis = UP - DOWN.
 2. Per perturbed gene g: effect_g = standardized (mean over g-cells - mean over NC cells)
    of the disease-axis score; same for module-24.
 3. Spearman rho between |effect_g| and atlas convergence_score (+10,000-perm null).
 4. Partial Spearman controlling for essentiality_chronos and log(n_cells).
 5. High-convergence (keep_in_pilot / top tier) vs low (drop_for_panel): Mann-Whitney U on |effect|.

Honest forward signal: screen not designed around our atlas.

Env: spatial (python/scanpy). Run via env binary directly.
"""
import os, json, time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import spearmanr, mannwhitneyu, rankdata

np.random.seed(42)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
H5AD = f"{ROOT}/data/perturbation/datasets/saunders2025/processed/saunders_perturb_hep_confident.h5ad"
DEG  = f"{ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"
MOD  = f"{ROOT}/Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv"
OV   = f"{ROOT}/data/perturbation/datasets/saunders2025/processed/saunders_panel_vs_masld_overlap.csv"
ORTH = f"{ROOT}/archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
ATLAS= f"{ROOT}/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
OUTDIR = f"{ROOT}/Analysis/SingleCell/results_gpu_v2/perturb_forward_validation"
os.makedirs(OUTDIR, exist_ok=True)
NC_LABEL = "control"

def log(*a): print(*a, flush=True)

# ---------- 1. Human signature + module-24 -> mouse ENSMUSG ----------
log("== Loading human signature ==")
deg = pd.read_csv(DEG)
deg["ensg"] = deg["gene"].str.replace(r"\.\d+$", "", regex=True)
sig = deg[(deg["padj"] < 0.05) & (deg["logFC"].abs() > 0.5)]
up_ensg = set(sig.loc[sig["logFC"] > 0, "ensg"])
dn_ensg = set(sig.loc[sig["logFC"] < 0, "ensg"])
log(f"disease-UP human={len(up_ensg)}  disease-DOWN human={len(dn_ensg)}")

mod = pd.read_csv(MOD, sep="\t")
m24_sym = set(mod.loc[mod["module"] == 24, "gene"])
# module genes are human SYMBOLS -> map to ENSG via DEG symbol table
sym2ensg = dict(zip(deg["symbol"], deg["ensg"]))
m24_ensg = {sym2ensg[s] for s in m24_sym if s in sym2ensg}
log(f"module-24 human symbols={len(m24_sym)}  mapped to ENSG={len(m24_ensg)}")

# human ENSG -> mouse ENSMUSG (direct Ensembl ortholog file: best coverage)
orth = pd.read_csv(ORTH, sep="\t").dropna(subset=["mouse_ensembl_gene_id", "human_ensembl_gene_id"])
h2m = orth.groupby("human_ensembl_gene_id")["mouse_ensembl_gene_id"].apply(list).to_dict()
def to_mouse(ensg_set):
    out = []
    for e in ensg_set:
        out.extend(h2m.get(e, []))
    return set(out)
up_mouse = to_mouse(up_ensg)
dn_mouse = to_mouse(dn_ensg)
m24_mouse = to_mouse(m24_ensg)
log(f"mouse-ortholog sets: UP={len(up_mouse)} DOWN={len(dn_mouse)} MOD24={len(m24_mouse)}")

# ---------- 2. Load mouse Perturb-seq (non-backed; backed mode hung) ----------
log("== Loading h5ad (non-backed) ==")
t0 = time.time()
A = ad.read_h5ad(H5AD)
log(f"loaded {A.shape} in {time.time()-t0:.0f}s; X dtype {A.X.dtype} max {A.X.max()}")
A.var_names = A.var_names.astype(str)  # ENSMUSG index
var_set = set(A.var_names)

def present(s):
    return [g for g in s if g in var_set]
up_in   = present(up_mouse)
dn_in   = present(dn_mouse)
m24_in  = present(m24_mouse)
log(f"in-h5ad signature genes: UP={len(up_in)} DOWN={len(dn_in)} MOD24={len(m24_in)}")

# normalize -> log1p (X is raw integer counts)
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)

# ---------- 3. Per-cell signature scoring (control-binned mean-z) ----------
log("== Scoring signatures (scanpy score_genes) ==")
sc.tl.score_genes(A, up_in,  score_name="score_up",   ctrl_size=50, random_state=42)
sc.tl.score_genes(A, dn_in,  score_name="score_dn",   ctrl_size=50, random_state=42)
sc.tl.score_genes(A, m24_in, score_name="score_mod24",ctrl_size=50, random_state=42)
A.obs["disease_axis"] = A.obs["score_up"] - A.obs["score_dn"]

obs = A.obs.copy()
obs["target_gene"] = obs["target_gene"].astype(str)
nc_mask = obs["target_gene"] == NC_LABEL
log(f"NC controls (label='{NC_LABEL}') n={int(nc_mask.sum())}")
nc_axis_mean = obs.loc[nc_mask, "disease_axis"].mean()
nc_axis_sd   = obs.loc[nc_mask, "disease_axis"].std(ddof=1)
nc_m24_mean  = obs.loc[nc_mask, "score_mod24"].mean()
nc_m24_sd    = obs.loc[nc_mask, "score_mod24"].std(ddof=1)
# pooled SD for standardization (use NC SD: stable, control-anchored)
log(f"NC disease_axis mean={nc_axis_mean:.4f} sd={nc_axis_sd:.4f}")

# ---------- 4. Per-perturbed-gene effect ----------
log("== Per-gene effects ==")
rows = []
for g, sub in obs.groupby("target_gene"):
    if g == NC_LABEL:
        continue
    n = len(sub)
    da_eff  = (sub["disease_axis"].mean()  - nc_axis_mean) / nc_axis_sd
    m24_eff = (sub["score_mod24"].mean()   - nc_m24_mean)  / nc_m24_sd
    hum = sub["human_ortholog"].astype(str)
    hum = hum[hum.notna() & (hum != "") & (hum != "nan")]
    human_ortholog = hum.mode().iloc[0] if len(hum) else np.nan
    rows.append(dict(target_gene=g, human_ortholog=human_ortholog, n_cells=n,
                     disease_axis_effect=da_eff, module24_effect=m24_eff))
res = pd.DataFrame(rows)
log(f"perturbed genes (excl NC): {len(res)}")

# ---------- 5. Join atlas convergence + essentiality ----------
ov = pd.read_csv(OV)
ov_keep = ov[["gene_mouse", "gene_human", "convergence_score", "convergence_rank",
              "tier", "recommendation", "coloc_best_susie_pp4", "is_conserved",
              "is_progression_driver"]].copy()
res = res.merge(ov_keep, left_on="target_gene", right_on="gene_mouse", how="left")
# fallback human ortholog from overlap CSV where obs was empty
res["human_ortholog"] = res["human_ortholog"].fillna(res["gene_human"])

atlas = pd.read_csv(ATLAS, usecols=["human_symbol", "essentiality_chronos"], low_memory=False)
atlas = atlas.dropna(subset=["human_symbol"]).drop_duplicates("human_symbol")
res = res.merge(atlas, left_on="human_ortholog", right_on="human_symbol", how="left")
res = res.drop(columns=["gene_mouse", "human_symbol"], errors="ignore")

res["abs_disease_axis_effect"] = res["disease_axis_effect"].abs()
res["abs_module24_effect"] = res["module24_effect"].abs()
res["log_n_cells"] = np.log(res["n_cells"])

# ---------- 6. Core statistics ----------
log("\n========== STATISTICS ==========")
core = res.dropna(subset=["convergence_score", "abs_disease_axis_effect"]).copy()
log(f"genes with convergence_score AND effect: {len(core)}")

def perm_spearman(x, y, n=10000, seed=42):
    rng = np.random.default_rng(seed)
    rho_obs = spearmanr(x, y).statistic
    yv = np.asarray(y)
    null = np.empty(n)
    for i in range(n):
        null[i] = spearmanr(x, rng.permutation(yv)).statistic
    p = (np.sum(np.abs(null) >= abs(rho_obs)) + 1) / (n + 1)
    return rho_obs, p

rho_raw, p_perm = perm_spearman(core["convergence_score"].values,
                                core["abs_disease_axis_effect"].values, n=10000)
log(f"RAW Spearman rho(|disease-axis effect|, convergence_score) = {rho_raw:.4f}")
log(f"  10,000-permutation empirical p = {p_perm:.4f}")

# module-24 raw
core24 = res.dropna(subset=["convergence_score", "abs_module24_effect"])
rho24, p24 = perm_spearman(core24["convergence_score"].values,
                           core24["abs_module24_effect"].values, n=10000)
log(f"RAW Spearman rho(|module-24 effect|, convergence_score) = {rho24:.4f}  perm p={p24:.4f}")

# ---------- 7. Partial Spearman controlling essentiality + log_n_cells ----------
def partial_spearman(df, x, y, covars):
    d = df.dropna(subset=[x, y] + covars).copy()
    if len(d) < 10:
        return np.nan, np.nan, len(d)
    # rank transform, regress out covars via OLS residuals, correlate residuals
    R = {c: rankdata(d[c].values) for c in [x, y] + covars}
    import numpy.linalg as la
    Z = np.column_stack([np.ones(len(d))] + [R[c] for c in covars])
    def resid(v):
        beta, *_ = la.lstsq(Z, v, rcond=None)
        return v - Z @ beta
    rx = resid(R[x]); ry = resid(R[y])
    # pearson on residual ranks = partial spearman
    rho = np.corrcoef(rx, ry)[0, 1]
    # permutation p on y residual labels
    rng = np.random.default_rng(7)
    null = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(10000)])
    p = (np.sum(np.abs(null) >= abs(rho)) + 1) / (10001)
    return rho, p, len(d)

rho_part, p_part, n_part = partial_spearman(
    res, "abs_disease_axis_effect", "convergence_score",
    ["essentiality_chronos", "log_n_cells"])
log(f"\nPARTIAL Spearman (control essentiality_chronos + log_n_cells):")
log(f"  rho = {rho_part:.4f}  perm p = {p_part:.4f}  (n={n_part})")

# also partial controlling only log_n_cells (power) and only essentiality
rho_pn, p_pn, _ = partial_spearman(res, "abs_disease_axis_effect", "convergence_score", ["log_n_cells"])
rho_pe, p_pe, _ = partial_spearman(res, "abs_disease_axis_effect", "convergence_score", ["essentiality_chronos"])
log(f"  partial (log_n_cells only):       rho={rho_pn:.4f} p={p_pn:.4f}")
log(f"  partial (essentiality only):      rho={rho_pe:.4f} p={p_pe:.4f}")

# ---------- 8. Stratified high vs low convergence ----------
log("\n== High vs low convergence (Mann-Whitney U on |disease-axis effect|) ==")
res["conv_group"] = np.where(
    res["recommendation"].isin(["keep_in_pilot", "keep_as_positive_control"]) |
    res["tier"].isin(["1_Genetic_validated", "2_Strong"]),
    "high", np.where(res["recommendation"] == "drop_for_panel", "low", "other"))
hi = res.loc[res["conv_group"] == "high", "abs_disease_axis_effect"].dropna()
lo = res.loc[res["conv_group"] == "low",  "abs_disease_axis_effect"].dropna()
log(f"high n={len(hi)} median|effect|={hi.median():.3f}  |  low n={len(lo)} median|effect|={lo.median():.3f}")
if len(hi) > 2 and len(lo) > 2:
    u, p_u = mannwhitneyu(hi, lo, alternative="greater")
    log(f"Mann-Whitney U (high > low, one-sided) U={u:.0f} p={p_u:.4f}")
    u2, p_u2 = mannwhitneyu(hi, lo, alternative="two-sided")
    log(f"Mann-Whitney U (two-sided) p={p_u2:.4f}")
else:
    p_u = np.nan; log("insufficient group sizes")

# ---------- 9. Save ----------
out_cols = ["target_gene", "human_ortholog", "convergence_score", "convergence_rank",
            "tier", "recommendation", "conv_group", "disease_axis_effect",
            "abs_disease_axis_effect", "module24_effect", "abs_module24_effect",
            "essentiality_chronos", "coloc_best_susie_pp4", "is_conserved",
            "is_progression_driver", "n_cells", "log_n_cells"]
res_out = res[[c for c in out_cols if c in res.columns]].sort_values(
    "convergence_score", ascending=False)
csv_path = f"{OUTDIR}/saunders_forward_validation_per_gene.csv"
res_out.to_csv(csv_path, index=False)
log(f"\nWrote per-gene CSV: {csv_path}")

stats = dict(
    n_perturbed_genes=int(len(res)),
    n_in_core=int(len(core)),
    n_NC_control_cells=int(nc_mask.sum()),
    n_total_cells=int(A.n_obs),
    sig_up_in_h5=len(up_in), sig_dn_in_h5=len(dn_in), mod24_in_h5=len(m24_in),
    rho_raw=float(rho_raw), perm_p_raw=float(p_perm),
    rho_module24=float(rho24), perm_p_module24=float(p24),
    rho_partial_ess_ncells=float(rho_part), perm_p_partial=float(p_part), n_partial=int(n_part),
    rho_partial_ncells_only=float(rho_pn), p_partial_ncells_only=float(p_pn),
    rho_partial_ess_only=float(rho_pe), p_partial_ess_only=float(p_pe),
    mwu_high_n=int(len(hi)), mwu_low_n=int(len(lo)),
    mwu_high_median=float(hi.median()), mwu_low_median=float(lo.median()),
    mwu_p_greater=float(p_u) if not np.isnan(p_u) else None,
)
with open(f"{OUTDIR}/forward_validation_stats.json", "w") as f:
    json.dump(stats, f, indent=2)
log(f"Wrote stats JSON: {OUTDIR}/forward_validation_stats.json")
log("\n==== DONE ====")
