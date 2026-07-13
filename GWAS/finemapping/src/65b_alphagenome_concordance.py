#!/usr/bin/env python3
"""
65b_alphagenome_concordance.py  --  AlphaGenome concordance reporting (Stage-2b).

Two independent read-outs, both additive (new file, script 65b):

  (A) AG vs measured Broadaway liver eQTL sign  -- the SAME pre-registered gate as
      64_concordance_benchmark.py (Borzoi), applied to AlphaGenome's signed gene
      LFC so the two models are judged on an identical yardstick:
        PRIMARY = sign-auROC of oriented ag_gene_lfc vs eqtl_sign on the
        is_signal_lead variants EXCLUDING strand_ambiguous palindromes,
        + 2000-boot 95% CI + 10,000-label-permutation null.
        TIERS: >=0.80&sig=STRONG ; 0.65-0.80&sig=MODEST ; <0.65 or not-sig=FAIL.

  (B) AG vs Borzoi two-model agreement -- the independent-second-model check. Both
      scorers run strictly hg38_ref -> hg38_alt, so raw sign(ag_gene_lfc) vs
      sign(borzoi_sed_rna) is directly comparable on rows whose hg38 alleles match.
      Reports % sign agreement + Spearman(rank) + a 2x2, on the common scored set.

ORIENTATION RULE (identical to 64):
  ag_gene_lfc is scored hg38_ref -> hg38_alt.
  effect_allele = truth ALT (eqtl_sign=+1 => effect allele raises expression).
    effect_allele == hg38_alt -> oriented =  ag_gene_lfc
    effect_allele == hg38_ref -> oriented = -ag_gene_lfc
    else (palindrome / mismatch) -> exclude.
"""
import os
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest
from sklearn.metrics import roc_auc_score

os.environ.setdefault("PYTHONHASHSEED", "0")
SEED = 42
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SEQ = f"{ROOT}/GWAS/finemapping/results/seqfunc"
AG = f"{SEQ}/alphagenome_eqtl_scores.tsv"
TRUTH = f"{SEQ}/broadaway_benchmark_truth.tsv"
BORZOI = f"{SEQ}/borzoi_eqtl_scores.tsv"
OUT_JSON = f"{SEQ}/AG_GATE_VERDICT.json"
OUT_PERROW = f"{SEQ}/ag_concordance_perrow.tsv"
OUT_AGVSBZ = f"{SEQ}/ag_vs_borzoi_agreement.json"
N_BOOT, N_PERM = 2000, 10000


def auroc(labels, scores):
    labels = np.asarray(labels); scores = np.asarray(scores)
    if len(np.unique(labels)) < 2:
        return np.nan
    return roc_auc_score(labels, scores)


def boot_ci_auroc(labels, scores, n=N_BOOT, seed=SEED):
    r = np.random.default_rng(seed)
    labels = np.asarray(labels); scores = np.asarray(scores); N = len(labels)
    vals = []; tries = 0
    while len(vals) < n and tries < n * 20:
        tries += 1
        idx = r.integers(0, N, N)
        if len(np.unique(labels[idx])) < 2:
            continue
        vals.append(roc_auc_score(labels[idx], scores[idx]))
    vals = np.array(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), vals


def perm_p_auroc(labels, scores, obs, n=N_PERM, seed=SEED):
    r = np.random.default_rng(seed + 1)
    labels = np.asarray(labels).copy(); scores = np.asarray(scores)
    perm = np.empty(n)
    for i in range(n):
        r.shuffle(labels)
        perm[i] = roc_auc_score(labels, scores)
    p = (1 + int(np.sum(perm >= obs))) / (1 + n)
    return float(p), perm


def strip_ver(e):
    return str(e).split(".")[0] if e is not None else e


# ---------------------------------------------------------------- load + join
ag = pd.read_csv(AG, sep="\t", dtype=str)
tr = pd.read_csv(TRUTH, sep="\t", dtype=str)

lead = tr[tr.is_signal_lead == "TRUE"].copy()
n_leads = len(lead)
n_palindrome = int((lead.strand_ambiguous == "TRUE").sum())
lead_np = lead[lead.strand_ambiguous == "FALSE"].copy()
n_after_palindrome = len(lead_np)

lead_np["key"] = lead_np.variant_id + "|" + lead_np.ensembl.map(strip_ver).fillna("")
agm = ag.copy()
agm["key"] = agm.variant_id + "|" + agm.ensembl.map(strip_ver).fillna("")
keep = ["key", "hg38_ref", "hg38_alt", "ag_gene_lfc", "ag_gene_quantile",
        "ag_splice_score", "ag_splice_quantile", "ag_accessibility_delta",
        "gene_matched", "is_indel", "note"]
m = lead_np.merge(agm[keep], on="key", how="left")
assert len(m) == n_after_palindrome, "join changed row count"

for c in ["ag_gene_lfc", "ag_gene_quantile", "ag_splice_score", "ag_accessibility_delta",
          "eqtl_beta", "eqtl_sign", "pp4_best", "eqtl_p"]:
    m[c] = pd.to_numeric(m[c], errors="coerce")


def orient(row, col):
    v = row[col]
    if pd.isna(v):
        return np.nan, "unscored"
    ea = row["effect_allele"]
    if ea == row["hg38_alt"]:
        return v, "direct"
    if ea == row["hg38_ref"]:
        return -v, "flip"
    return np.nan, "mismatch"


oc = m.apply(lambda r: orient(r, "ag_gene_lfc"), axis=1)
m["oriented_lfc"] = [t[0] for t in oc]
m["orient_case"] = [t[1] for t in oc]
occ = m.apply(lambda r: orient(r, "ag_accessibility_delta"), axis=1)
m["oriented_acc"] = [t[0] for t in occ]
m["label"] = (m["eqtl_sign"] > 0).astype(int)

b = m[m["oriented_lfc"].notna()].copy()
n_bench = len(b)
n_unscored = n_after_palindrome - n_bench
print(f"leads={n_leads} palindromes_excluded={n_palindrome} "
      f"after_palindrome={n_after_palindrome} unscored_excluded={n_unscored} n_benchmark={n_bench}")
print("orient_case:", b.orient_case.value_counts().to_dict())
print("label balance (1=+1 / 0=-1):", b.label.value_counts().to_dict())

# ---------------------------------------------------------------- (A) PRIMARY sign-auROC
labels = b["label"].values
scores = b["oriented_lfc"].values
auc = auroc(labels, scores)
ci_lo, ci_hi, _ = boot_ci_auroc(labels, scores)
perm_p, perm_vals = perm_p_auroc(labels, scores, auc)
concord = float(np.mean(np.sign(scores) == np.where(labels == 1, 1, -1)) * 100)
n_conc = int(round(concord / 100 * n_bench))
conc_binom_p = binomtest(n_conc, n_bench, 0.5, alternative="greater").pvalue
sp_signed_rho, sp_signed_p = spearmanr(b["oriented_lfc"].values, b["eqtl_beta"].values)

# accessibility secondary
acc = b[b["oriented_acc"].notna()].copy()
if len(acc) and len(np.unique(acc["label"])) == 2:
    acc_auc = auroc(acc["label"].values, acc["oriented_acc"].values)
else:
    acc_auc = np.nan

sig = perm_p < 0.05
if np.isnan(auc):
    tier = "FAIL"
elif auc >= 0.80 and sig:
    tier = "STRONG"
elif auc >= 0.65 and sig:
    tier = "MODEST"
else:
    tier = "FAIL"

# anchor QC
anc = ag[ag.is_anchor == "TRUE"]
sort1 = anc[anc.anchor_name == "SORT1_posctrl"]
sort1_lfc = float(sort1.ag_gene_lfc.iloc[0]) if len(sort1) else np.nan
hsd = anc[anc.anchor_name == "HSD17B13_splice"]
hsd_splice_q = float(hsd.ag_splice_quantile.iloc[0]) if len(hsd) else np.nan

gate = dict(
    model="AlphaGenome (GeneMaskLFCScorer, liver RNA_SEQ, DIFF_LOG2_SUM)",
    auc=round(float(auc), 4), ci=[round(ci_lo, 4), round(ci_hi, 4)],
    ci_method="2000-bootstrap percentile 95%",
    p=float(perm_p), p_method="10000-label-permutation, one-sided (auROC>=obs)",
    significant=bool(sig), tier=tier, n=int(n_bench),
    sign_concordance_pct=round(concord, 2), concordance_binom_p=float(conc_binom_p),
    spearman_signed_oriented_lfc_vs_beta=round(float(sp_signed_rho), 4),
    accessibility_auc=(None if np.isnan(acc_auc) else round(float(acc_auc), 4)),
    accounting=dict(leads_total=n_leads, palindromes_excluded=n_palindrome,
                    unscored_excluded=int(n_unscored), n_benchmark=int(n_bench)),
    anchor_SORT1_ag_gene_lfc=(None if np.isnan(sort1_lfc) else round(sort1_lfc, 5)),
    anchor_SORT1_pass=(bool(sort1_lfc > 0) if not np.isnan(sort1_lfc) else None),
    anchor_HSD17B13_splice_quantile=(None if np.isnan(hsd_splice_q) else round(hsd_splice_q, 5)),
    label_balance=dict(pos_plus1=int((b.label == 1).sum()), neg_minus1=int((b.label == 0).sum())),
    seed=SEED,
    tier_rule=">=0.80&sig=STRONG ; 0.65-0.80&sig=MODEST ; <0.65 or not-sig=FAIL",
)
with open(OUT_JSON, "w") as f:
    json.dump(gate, f, indent=2)

b_out = b[["variant_id", "gene", "ensembl", "effect_allele", "hg38_ref", "hg38_alt",
           "orient_case", "ag_gene_lfc", "oriented_lfc", "eqtl_beta", "eqtl_sign",
           "label", "pp4_best", "study", "trait", "ancestry"]].copy()
b_out["lfc_sign"] = np.sign(b_out["oriented_lfc"])
b_out["concordant"] = (b_out["lfc_sign"] == np.where(b_out["label"] == 1, 1, -1))
b_out.to_csv(OUT_PERROW, sep="\t", index=False)

# ---------------------------------------------------------------- (B) AG vs Borzoi
agbz = dict(available=False)
if os.path.exists(BORZOI):
    bz = pd.read_csv(BORZOI, sep="\t", dtype=str)
    bz["key"] = bz.variant_id + "|" + bz.ensembl.map(strip_ver).fillna("")
    agm2 = ag.copy()
    agm2["key"] = agm2.variant_id + "|" + agm2.ensembl.map(strip_ver).fillna("")
    j = agm2.merge(bz[["key", "hg38_ref", "hg38_alt", "borzoi_sed_rna"]],
                   on="key", how="inner", suffixes=("_ag", "_bz"))
    j["ag_gene_lfc"] = pd.to_numeric(j["ag_gene_lfc"], errors="coerce")
    j["borzoi_sed_rna"] = pd.to_numeric(j["borzoi_sed_rna"], errors="coerce")
    # both scored hg38 ref->alt: require identical hg38 alleles for a clean sign compare
    allele_ok = (j["hg38_ref_ag"].astype(str) == j["hg38_ref_bz"].astype(str)) & \
                (j["hg38_alt_ag"].astype(str) == j["hg38_alt_bz"].astype(str))
    c = j[allele_ok & j.ag_gene_lfc.notna() & j.borzoi_sed_rna.notna()].copy()
    if len(c):
        sa = float(np.mean(np.sign(c.ag_gene_lfc) == np.sign(c.borzoi_sed_rna)) * 100)
        n_agree = int(round(sa / 100 * len(c)))
        agree_p = binomtest(n_agree, len(c), 0.5, alternative="greater").pvalue
        rho, rho_p = spearmanr(c.ag_gene_lfc, c.borzoi_sed_rna)
        both_pos = int(((c.ag_gene_lfc > 0) & (c.borzoi_sed_rna > 0)).sum())
        both_neg = int(((c.ag_gene_lfc < 0) & (c.borzoi_sed_rna < 0)).sum())
        ag_pos_bz_neg = int(((c.ag_gene_lfc > 0) & (c.borzoi_sed_rna < 0)).sum())
        ag_neg_bz_pos = int(((c.ag_gene_lfc < 0) & (c.borzoi_sed_rna > 0)).sum())
        agbz = dict(available=True, n_common_scored=int(len(c)),
                    n_allele_mismatch_excluded=int((~allele_ok).sum()),
                    sign_agreement_pct=round(sa, 2), agreement_binom_p=float(agree_p),
                    spearman_rho=round(float(rho), 4), spearman_p=float(rho_p),
                    contingency=dict(both_pos=both_pos, both_neg=both_neg,
                                     ag_pos_borzoi_neg=ag_pos_bz_neg,
                                     ag_neg_borzoi_pos=ag_neg_bz_pos),
                    note="raw hg38 ref->alt sign agreement (pre-orientation); "
                         "magnitudes differ by scale (AG log2-LFC vs Borzoi SED) so rank/sign used")
with open(OUT_AGVSBZ, "w") as f:
    json.dump(agbz, f, indent=2)

# ---------------------------------------------------------------- (C) matched head-to-head
# AG vs Borzoi sign-auROC on the IDENTICAL benchmark leads (the ~186-lead set both
# models scored). Orient Borzoi's SED with the same effect-allele contract, intersect
# with AG's oriented set, and score BOTH models on exactly those variants -- the fair
# apples-to-apples head-to-head requested by the team lead.
OUT_MATCHED = f"{SEQ}/ag_vs_borzoi_matched_auroc.json"
matched = dict(available=False)
if os.path.exists(BORZOI):
    bz2 = pd.read_csv(BORZOI, sep="\t", dtype=str)
    bz2["key"] = bz2.variant_id + "|" + bz2.ensembl.map(strip_ver).fillna("")
    bz2 = bz2.rename(columns={"hg38_ref": "bz_hg38_ref", "hg38_alt": "bz_hg38_alt"})
    bz2["borzoi_sed_rna"] = pd.to_numeric(bz2["borzoi_sed_rna"], errors="coerce")
    mm = m.merge(bz2[["key", "bz_hg38_ref", "bz_hg38_alt", "borzoi_sed_rna"]],
                 on="key", how="left")

    def orient_bz(row):
        sed = row["borzoi_sed_rna"]
        if pd.isna(sed):
            return np.nan
        ea = row["effect_allele"]
        if ea == row["bz_hg38_alt"]:
            return sed
        if ea == row["bz_hg38_ref"]:
            return -sed
        return np.nan

    mm["oriented_bz"] = mm.apply(orient_bz, axis=1)
    mset = mm[mm["oriented_lfc"].notna() & mm["oriented_bz"].notna()].copy()
    if len(mset) and mset["label"].nunique() == 2:
        lab = mset["label"].values
        ag_sc = mset["oriented_lfc"].values
        bz_sc = mset["oriented_bz"].values
        ag_auc = auroc(lab, ag_sc); bz_auc = auroc(lab, bz_sc)
        ag_lo, ag_hi, _ = boot_ci_auroc(lab, ag_sc)
        bz_lo, bz_hi, _ = boot_ci_auroc(lab, bz_sc)
        ag_pp, _ = perm_p_auroc(lab, ag_sc, ag_auc)
        bz_pp, _ = perm_p_auroc(lab, bz_sc, bz_auc)
        ag_conc = float(np.mean(np.sign(ag_sc) == np.where(lab == 1, 1, -1)) * 100)
        bz_conc = float(np.mean(np.sign(bz_sc) == np.where(lab == 1, 1, -1)) * 100)

        def tier_of(a, p):
            if np.isnan(a):
                return "FAIL"
            if a >= 0.80 and p < 0.05:
                return "STRONG"
            if a >= 0.65 and p < 0.05:
                return "MODEST"
            return "FAIL"

        matched = dict(
            available=True, n_matched=int(len(mset)),
            label_balance=dict(pos=int((lab == 1).sum()), neg=int((lab == 0).sum())),
            alphagenome=dict(auroc=round(float(ag_auc), 4), ci=[round(ag_lo, 4), round(ag_hi, 4)],
                             perm_p=float(ag_pp), concordance_pct=round(ag_conc, 2),
                             tier=tier_of(ag_auc, ag_pp)),
            borzoi=dict(auroc=round(float(bz_auc), 4), ci=[round(bz_lo, 4), round(bz_hi, 4)],
                        perm_p=float(bz_pp), concordance_pct=round(bz_conc, 2),
                        tier=tier_of(bz_auc, bz_pp)),
            note="identical benchmark leads; both models oriented to effect allele; same "
                 "pre-registered tiers (>=0.80 STRONG / 0.65-0.80&sig MODEST / else FAIL)")
with open(OUT_MATCHED, "w") as f:
    json.dump(matched, f, indent=2)

# ---------------------------------------------------------------- console
print("\n============== (A) AG vs Broadaway eQTL ==============")
print(f"n_benchmark        : {n_bench} (+1:{int((b.label==1).sum())} / -1:{int((b.label==0).sum())})")
print(f"PRIMARY sign-auROC : {auc:.4f}  95% CI [{ci_lo:.4f}, {ci_hi:.4f}]")
print(f"permutation p      : {perm_p:.5f}  (sig={sig})")
print(f"sign-concordance   : {concord:.1f}%  (binom p {conc_binom_p:.3g})")
print(f"spearman signed    : {sp_signed_rho:.3f}  p={sp_signed_p:.3g}")
print(f"accessibility auROC: {acc_auc if np.isnan(acc_auc) else round(acc_auc,4)}")
print(f"SORT1 anchor lfc   : {sort1_lfc:+.5f}  (must be >0)")
print(f"==> TIER           : {tier}")
print("\n============== (B) AG vs Borzoi (two-model) ==========")
if agbz.get("available"):
    print(f"n common scored    : {agbz['n_common_scored']}")
    print(f"sign agreement     : {agbz['sign_agreement_pct']}%  (binom p {agbz['agreement_binom_p']:.3g})")
    print(f"spearman rho       : {agbz['spearman_rho']}  p={agbz['spearman_p']:.3g}")
    print(f"contingency        : {agbz['contingency']}")
else:
    print("borzoi_eqtl_scores.tsv not available -> AG-vs-Borzoi skipped")
print("\n===== (C) MATCHED head-to-head (identical benchmark leads) =====")
if matched.get("available"):
    a, bb = matched["alphagenome"], matched["borzoi"]
    print(f"n_matched          : {matched['n_matched']}  "
          f"(+1:{matched['label_balance']['pos']} / -1:{matched['label_balance']['neg']})")
    print(f"AlphaGenome auROC  : {a['auroc']:.4f} {a['ci']}  perm p={a['perm_p']:.4f}  "
          f"conc {a['concordance_pct']}%  TIER {a['tier']}")
    print(f"Borzoi     auROC   : {bb['auroc']:.4f} {bb['ci']}  perm p={bb['perm_p']:.4f}  "
          f"conc {bb['concordance_pct']}%  TIER {bb['tier']}")
else:
    print("matched head-to-head unavailable")
print("\nwrote:", OUT_JSON, "|", OUT_PERROW, "|", OUT_AGVSBZ, "|", OUT_MATCHED)
