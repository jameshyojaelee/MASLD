#!/usr/bin/env python
"""
64_concordance_benchmark.py  —  Borzoi vs measured Broadaway liver eQTL sign benchmark.

Pre-registered gate (locked before any number was seen):
  PRIMARY = sign-auROC of Borzoi signed gene-level SED (oriented to the eQTL EFFECT allele)
  vs measured eQTL sign, on the 234 is_signal_lead variants EXCLUDING strand_ambiguous
  palindromes, with 2000-bootstrap 95% CI + 10,000-label-permutation null.
  TIERS:  >=0.80 & sig = STRONG ; 0.65-0.80 & sig = MODEST (nomination-only) ;
          <0.65 or not-sig = FAIL.

ORIENTATION RULE (the sign-flip contract), applied EXACTLY:
  Borzoi borzoi_sed_rna is scored hg38_ref -> hg38_alt.
  eqtl_sign is oriented to effect_allele = truth ALT (eqtl_sign=+1 => ALT raises expr).
    if effect_allele == hg38_alt  -> oriented_sed =  borzoi_sed_rna   (direct)
    if effect_allele == hg38_ref  -> oriented_sed = -borzoi_sed_rna   (flip)
    else (palindrome / unscored) -> exclude.
  Verified on SORT1 anchor: effect_allele T == hg38_alt T => direct; SED +28.39 & sign +1 => concordant.

ADDITIVE ONLY: reads Phase-1 outputs, writes results/seqfunc/{concordance_benchmark_results.csv,
GATE_VERDICT.json, concordance_benchmark_perrow.tsv}. Does not edit any existing script.
"""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

os.environ.setdefault("PYTHONHASHSEED", "0")
SEED = 42
rng = np.random.default_rng(SEED)

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SEQ  = f"{ROOT}/GWAS/finemapping/results/seqfunc"
SCORES = f"{SEQ}/borzoi_eqtl_scores.tsv"
TRUTH  = f"{SEQ}/broadaway_benchmark_truth.tsv"
OUT_CSV   = f"{SEQ}/concordance_benchmark_results.csv"
OUT_JSON  = f"{SEQ}/GATE_VERDICT.json"
OUT_PERROW= f"{SEQ}/concordance_benchmark_perrow.tsv"

N_BOOT = 2000
N_PERM = 10000

# ---------------------------------------------------------------- helpers
def auroc(labels, scores):
    """sign-auROC: labels in {0,1} (1 = eqtl_sign +1), scores = oriented SED."""
    labels = np.asarray(labels); scores = np.asarray(scores)
    if len(np.unique(labels)) < 2:
        return np.nan
    return roc_auc_score(labels, scores)

def boot_ci_auroc(labels, scores, n=N_BOOT, seed=SEED):
    r = np.random.default_rng(seed)
    labels = np.asarray(labels); scores = np.asarray(scores); N = len(labels)
    vals = []
    tries = 0
    while len(vals) < n and tries < n * 20:
        tries += 1
        idx = r.integers(0, N, N)
        if len(np.unique(labels[idx])) < 2:
            continue
        vals.append(roc_auc_score(labels[idx], scores[idx]))
    vals = np.array(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), vals

def perm_p_auroc(labels, scores, obs, n=N_PERM, seed=SEED):
    """One-sided label-permutation p: P(perm auROC >= observed) under H0 (no info)."""
    r = np.random.default_rng(seed + 1)
    labels = np.asarray(labels).copy(); scores = np.asarray(scores)
    perm = np.empty(n)
    for i in range(n):
        r.shuffle(labels)
        perm[i] = roc_auc_score(labels, scores)
    p = (1 + int(np.sum(perm >= obs))) / (1 + n)
    return float(p), perm

def boot_ci_spearman(x, y, n=N_BOOT, seed=SEED):
    r = np.random.default_rng(seed + 2)
    x = np.asarray(x); y = np.asarray(y); N = len(x)
    vals = []
    for _ in range(n):
        idx = r.integers(0, N, N)
        rho, _ = spearmanr(x[idx], y[idx])
        if not np.isnan(rho):
            vals.append(rho)
    vals = np.array(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))

# ---------------------------------------------------------------- load + join
sc = pd.read_csv(SCORES, sep="\t", dtype=str)
tr = pd.read_csv(TRUTH,  sep="\t", dtype=str)

lead = tr[tr.is_signal_lead == "TRUE"].copy()
n_leads = len(lead)
n_palindrome = int((lead.strand_ambiguous == "TRUE").sum())
lead_np = lead[lead.strand_ambiguous == "FALSE"].copy()      # exclude palindromes
n_after_palindrome = len(lead_np)

lead_np["key"] = lead_np.variant_id + "|" + lead_np.ensembl.fillna("")
scm = sc.copy()
scm["key"] = scm.variant_id + "|" + scm.ensembl.fillna("")
scm_keep = scm[["key", "hg38_ref", "hg38_alt", "borzoi_sed_rna",
                "borzoi_sed_rna_alltracks", "borzoi_acc_delta",
                "gene_in_window", "is_indel", "n_rna_tracks_used", "note"]]
m = lead_np.merge(scm_keep, on="key", how="left")
assert len(m) == n_after_palindrome, "join changed row count"

# numeric casts
for c in ["borzoi_sed_rna", "borzoi_sed_rna_alltracks", "borzoi_acc_delta",
          "eqtl_beta", "eqtl_sign", "pp4_best", "eqtl_p"]:
    m[c] = pd.to_numeric(m[c], errors="coerce")

# ---------------------------------------------------------------- orientation
def orient(row, sedcol):
    sed = row[sedcol]
    if pd.isna(sed):
        return np.nan, "unscored"
    ea = row["effect_allele"]
    if ea == row["hg38_alt"]:
        return sed, "direct"
    if ea == row["hg38_ref"]:
        return -sed, "flip"
    return np.nan, "mismatch"

oc = m.apply(lambda r: orient(r, "borzoi_sed_rna"), axis=1)
m["oriented_sed"]  = [t[0] for t in oc]
m["orient_case"]   = [t[1] for t in oc]
oca = m.apply(lambda r: orient(r, "borzoi_sed_rna_alltracks"), axis=1)
m["oriented_sed_alltracks"] = [t[0] for t in oca]
occ = m.apply(lambda r: orient(r, "borzoi_acc_delta"), axis=1)
m["oriented_acc"] = [t[0] for t in occ]

# binary label: 1 = eqtl_sign +1 (ALT/effect allele increases expression)
m["label"] = (m["eqtl_sign"] > 0).astype(int)

# in-window scored set = usable oriented SED
b = m[m["oriented_sed"].notna()].copy()
n_bench = len(b)
n_unscored = n_after_palindrome - n_bench

print(f"leads={n_leads} palindromes_excluded={n_palindrome} "
      f"after_palindrome={n_after_palindrome} unscored_excluded={n_unscored} n_benchmark={n_bench}")
print("orient_case (benchmark set):", b.orient_case.value_counts().to_dict())
print("label balance (1=+1 / 0=-1):", b.label.value_counts().to_dict())

# ---------------------------------------------------------------- PRIMARY sign-auROC
labels = b["label"].values
scores = b["oriented_sed"].values
auc = auroc(labels, scores)
ci_lo, ci_hi, boot_vals = boot_ci_auroc(labels, scores)
perm_p, perm_vals = perm_p_auroc(labels, scores, auc)
concord = float(np.mean(np.sign(scores) == np.where(labels == 1, 1, -1)) * 100)
# binomial p for concordance vs 0.5 (two-sided, normal approx)
from scipy.stats import binomtest
n_conc = int(round(concord / 100 * n_bench))
conc_binom_p = binomtest(n_conc, n_bench, 0.5, alternative="greater").pvalue

# alltracks secondary auROC
ba = b[b["oriented_sed_alltracks"].notna()]
auc_alltracks = auroc(ba["label"].values, ba["oriented_sed_alltracks"].values)

# ---------------------------------------------------------------- MAGNITUDE Spearman
# signed: oriented SED vs eqtl_beta (both oriented to effect allele)
sp_signed_rho, sp_signed_p = spearmanr(b["oriented_sed"].values, b["eqtl_beta"].values)
sp_signed_lo, sp_signed_hi = boot_ci_spearman(b["oriented_sed"].values, b["eqtl_beta"].values)
# magnitude: |SED| vs |beta|
absSED = np.abs(b["borzoi_sed_rna"].values); absB = np.abs(b["eqtl_beta"].values)
sp_mag_rho, sp_mag_p = spearmanr(absSED, absB)
sp_mag_lo, sp_mag_hi = boot_ci_spearman(absSED, absB)

# ---------------------------------------------------------------- STRATIFIED auROC
def tertile_auroc(df, col):
    q = df[col].quantile([1/3, 2/3]).values
    out = {}
    for name, mask in [("low", df[col] <= q[0]),
                       ("mid", (df[col] > q[0]) & (df[col] <= q[1])),
                       ("high", df[col] > q[1])]:
        sub = df[mask]
        out[name] = (len(sub), auroc(sub["label"].values, sub["oriented_sed"].values))
    return q, out

pp4_q, pp4_strat = tertile_auroc(b, "pp4_best")
b["_absbeta"] = np.abs(b["eqtl_beta"])
absb_q, absb_strat = tertile_auroc(b, "_absbeta")

# ---------------------------------------------------------------- ACCESSIBILITY (secondary)
acc = b[b["oriented_acc"].notna()].copy()
n_acc = len(acc)
if n_acc > 0 and len(np.unique(acc["label"])) == 2:
    acc_auc = auroc(acc["label"].values, acc["oriented_acc"].values)
    acc_conc = float(np.mean(np.sign(acc["oriented_acc"].values) ==
                             np.where(acc["label"].values == 1, 1, -1)) * 100)
    acc_lo, acc_hi, _ = boot_ci_auroc(acc["label"].values, acc["oriented_acc"].values)
    acc_perm_p, _ = perm_p_auroc(acc["label"].values, acc["oriented_acc"].values, acc_auc)
else:
    acc_auc = acc_conc = acc_lo = acc_hi = acc_perm_p = np.nan

# ---------------------------------------------------------------- SORT1 anchor QC
anc = sc[sc.is_anchor == "TRUE"].iloc[0]
anc_sed = float(anc.borzoi_sed_rna)

# ---------------------------------------------------------------- TIER decision
sig = perm_p < 0.05
if np.isnan(auc):
    tier = "FAIL"
elif auc >= 0.80 and sig:
    tier = "STRONG"
elif auc >= 0.65 and sig:
    tier = "MODEST"
else:
    tier = "FAIL"

# ---------------------------------------------------------------- write results CSV
rows = []
def add(metric, stratum, n, est, lo, hi, p, note=""):
    rows.append(dict(metric=metric, stratum=stratum, n=n,
                     estimate=est, ci_low=lo, ci_high=hi, p_value=p, note=note))

add("sign_auroc", "primary_liver_SED", n_bench, auc, ci_lo, ci_hi, perm_p,
    "oriented Borzoi liver-RNA SED vs eqtl_sign; 2000-boot CI; 10000-perm p (one-sided)")
add("sign_concordance_pct", "primary_liver_SED", n_bench, concord, np.nan, np.nan, conc_binom_p,
    "sign(oriented SED)==eqtl_sign; binomial p vs 0.5 (greater)")
add("sign_auroc", "secondary_alltracks_SED", len(ba), auc_alltracks, np.nan, np.nan, np.nan,
    "oriented Borzoi all-tracks SED (robustness)")
add("spearman_signed", "oriented_SED_vs_eqtl_beta", n_bench, sp_signed_rho, sp_signed_lo, sp_signed_hi, sp_signed_p,
    "signed magnitude concordance; 2000-boot CI")
add("spearman_magnitude", "absSED_vs_absbeta", n_bench, sp_mag_rho, sp_mag_lo, sp_mag_hi, sp_mag_p,
    "|SED| vs |beta|; 2000-boot CI")
for name in ["low", "mid", "high"]:
    n_s, a_s = pp4_strat[name]
    add("sign_auroc", f"pp4_tertile_{name}", n_s, a_s, np.nan, np.nan, np.nan,
        f"pp4_best tertile cutpoints={np.round(pp4_q,4).tolist()}")
for name in ["low", "mid", "high"]:
    n_s, a_s = absb_strat[name]
    add("sign_auroc", f"absbeta_tertile_{name}", n_s, a_s, np.nan, np.nan, np.nan,
        f"|eqtl_beta| tertile cutpoints={np.round(absb_q,4).tolist()}")
add("sign_auroc", "accessibility_delta", n_acc, acc_auc, acc_lo, acc_hi, acc_perm_p,
    "oriented Borzoi accessibility SED vs eqtl_sign (secondary)")
add("sign_concordance_pct", "accessibility_delta", n_acc, acc_conc, np.nan, np.nan, np.nan,
    "sign(oriented acc)==eqtl_sign")
add("qc_anchor_SORT1_SED", "positive_control", 1, anc_sed, np.nan, np.nan, np.nan,
    "SORT1 rs12740374 G->T must be POSITIVE (Musunuru 2010)")
add("n_accounting", "leads_total", n_leads, np.nan, np.nan, np.nan, np.nan, "234 independent signal leads")
add("n_accounting", "palindromes_excluded", n_palindrome, np.nan, np.nan, np.nan, np.nan, "strand_ambiguous==TRUE")
add("n_accounting", "unscored_excluded", n_unscored, np.nan, np.nan, np.nan, np.nan, "no Borzoi SED / outside receptive field")
add("n_accounting", "n_benchmark", n_bench, np.nan, np.nan, np.nan, np.nan, "scored, oriented, non-palindrome")

res = pd.DataFrame(rows)
res.to_csv(OUT_CSV, index=False)

# per-row provenance
b_out = b[["variant_id", "gene", "ensembl", "effect_allele", "hg38_ref", "hg38_alt",
           "orient_case", "borzoi_sed_rna", "oriented_sed", "eqtl_beta", "eqtl_sign",
           "label", "pp4_best", "eqtl_p", "study", "trait", "ancestry", "oriented_acc"]].copy()
b_out["sed_sign"] = np.sign(b_out["oriented_sed"])
b_out["concordant"] = (b_out["sed_sign"] == np.where(b_out["label"] == 1, 1, -1))
b_out.to_csv(OUT_PERROW, sep="\t", index=False)

# ---------------------------------------------------------------- GATE json
gate = dict(
    auc=round(float(auc), 4),
    ci=[round(ci_lo, 4), round(ci_hi, 4)],
    ci_method="2000-bootstrap percentile 95%",
    p=float(perm_p),
    p_method="10000-label-permutation, one-sided (auROC>=obs)",
    perm_null_median=round(float(np.median(perm_vals)), 4),
    significant=bool(sig),
    tier=tier,
    n=int(n_bench),
    sign_concordance_pct=round(concord, 2),
    concordance_binom_p=float(conc_binom_p),
    spearman_signed=round(float(sp_signed_rho), 4),
    spearman_signed_ci=[round(sp_signed_lo, 4), round(sp_signed_hi, 4)],
    spearman_magnitude=round(float(sp_mag_rho), 4),
    auc_alltracks=round(float(auc_alltracks), 4),
    accessibility_auc=(round(float(acc_auc), 4) if not np.isnan(acc_auc) else None),
    accessibility_n=int(n_acc),
    accounting=dict(leads_total=n_leads, palindromes_excluded=n_palindrome,
                    unscored_excluded=int(n_unscored), n_benchmark=int(n_bench)),
    anchor_SORT1_SED=round(anc_sed, 4),
    label_balance=dict(pos_plus1=int((b.label == 1).sum()), neg_minus1=int((b.label == 0).sum())),
    seed=SEED,
    tier_rule=">=0.80&sig=STRONG ; 0.65-0.80&sig=MODEST ; <0.65 or not-sig=FAIL",
)
with open(OUT_JSON, "w") as f:
    json.dump(gate, f, indent=2)

# ---------------------------------------------------------------- console summary
print("\n================ GATE ================")
print(f"n_benchmark            : {n_bench}  (label +1: {int((b.label==1).sum())} / -1: {int((b.label==0).sum())})")
print(f"PRIMARY sign-auROC     : {auc:.4f}  95% CI [{ci_lo:.4f}, {ci_hi:.4f}]")
print(f"permutation p (1-sided): {perm_p:.5f}  (null median {np.median(perm_vals):.3f})")
print(f"sign-concordance       : {concord:.1f}%  (binom p {conc_binom_p:.3g})")
print(f"spearman signed        : {sp_signed_rho:.3f} [{sp_signed_lo:.3f},{sp_signed_hi:.3f}]  p={sp_signed_p:.3g}")
print(f"spearman |SED| vs |b|  : {sp_mag_rho:.3f} [{sp_mag_lo:.3f},{sp_mag_hi:.3f}]  p={sp_mag_p:.3g}")
print(f"auROC alltracks        : {auc_alltracks:.4f}")
print(f"pp4 tertiles (n,auc)   : " + " ".join(f"{k}={pp4_strat[k][0]},{pp4_strat[k][1]:.3f}" for k in ['low','mid','high']))
print(f"|beta| tertiles (n,auc): " + " ".join(f"{k}={absb_strat[k][0]},{absb_strat[k][1]:.3f}" for k in ['low','mid','high']))
print(f"accessibility auROC    : {acc_auc if np.isnan(acc_auc) else round(acc_auc,4)}  (n={n_acc}, concord {acc_conc if np.isnan(acc_conc) else round(acc_conc,1)}%)")
print(f"SORT1 anchor SED       : {anc_sed:.3f}  (must be >0)")
print(f"==> TIER               : {tier}")
print("======================================")
print(f"wrote {OUT_CSV}")
print(f"wrote {OUT_JSON}")
print(f"wrote {OUT_PERROW}")
