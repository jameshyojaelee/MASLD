#!/usr/bin/env python3
"""
85_eqtl_signal_diagnose.py — Phase-6 WS-2: signal-level, ToS-clean eQTL DIRECTION
diagnostic (supersedes the eQTL arm of 80_direction_hardened.py).

WHY (red-team critique F3/M1, source-verified):
  * The old eQTL panel scored 211 "independent leads" that collapse to only
    ~182 unique physical positions / ~169 independent (variant,gene) signals —
    one position appears 5x. Row-bootstrap over these UNDERSTATES uncertainty.
  * The old headline trained an elastic-net on AlphaGenome output features (ToS
    violation) and framed the result as "SOTA models FAIL eQTL direction."

FIX (this script):
  * ZERO-SHOT ONLY, OPEN license: Borzoi logSED (johahi/borzoi-pytorch, 64b) is
    the primary. NO trained model, NO AlphaGenome feature used here (the AG-native
    zero-shot diagnostic is a SEPARATE run, 85b, reported with the ToS notice).
  * Signal-level units: collapse to one row per Broadaway signal_id (the (variant,
    gene) eQTL signal) AND report a one-per-physical-variant sensitivity. CIs are
    block-bootstrapped by whole CHROMOSOME (LD-robust); a signal-clustered
    bootstrap is also reported.
  * strand_ambiguous already dropped upstream (verified 0/211); we re-assert it.
  * FRAMING: this is a diagnostic, not a verdict on the models. Report as
    "in this Broadaway/colocalization-derived liver benchmark, zero-shot Borzoi
    logSED sign-auROC did NOT meet the prespecified 0.65 threshold" — NOT
    "the models fail eQTL direction." The point of 85b (AG-native) is to DIAGNOSE
    whether this is ascertainment/dependence vs a genuinely harder liver panel
    (cf. AlphaGenome's published ~0.80 on fine-mapped GTEx).

APPLY-ONLY FIREWALL (unchanged): benchmark only; never a per-locus direction call
  on eQTL-absent/convergence loci; never summed into convergence.
"""
import json
import os
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
LAB = os.path.join(SF, "direction_labels_eqtl.tsv")
TRUTH = os.path.join(SF, "broadaway_benchmark_truth.tsv")
BZ = os.path.join(SF, "borzoi_eqtl_logsed_scores.tsv")
OUT_JSON = os.path.join(SF, "eqtl_signal_diagnose_verdict.json")
OUT_TSV = os.path.join(SF, "eqtl_signal_diagnose_permodel.tsv")

PRIMARY_BAR = 0.65
NBOOT = 2000
NPERM = 1000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
ZS_MODELS = {                       # zero-shot, open-license signed deltas
    "borzoi_logsed_liver": "borzoi_logsed_liver",
    "borzoi_logsed_allliver": "borzoi_logsed_allliver",
    "borzoi_acc_delta": "borzoi_acc_delta",
}
AG_ZS = {                           # AlphaGenome native eQTL scorer (GeneMaskLFC), ZERO-SHOT ANNOTATION ONLY
    "alphagenome_gene_lfc": "ag_gene_lfc",
    "alphagenome_gene_lfc_alltracks": "ag_gene_lfc_alltracks",
}
AG_FILE = os.path.join(SF, "alphagenome_eqtl_scores.tsv")
AG_NOTICE = ("AlphaGenome ag_gene_lfc is its NATIVE eQTL scorer (GeneMaskLFCScorer); used here "
             "strictly as ZERO-SHOT annotation/evaluation, NOT to train any model (AlphaGenome FAQ "
             "output restriction; non-commercial). Diagnostic purpose: if AG is ALSO ~chance on this "
             "panel, the sub-threshold result reflects ascertainment (independent liver-microarray "
             "marginal leads, out of the fine-mapped-GTEx distribution where AG reports ~0.80), not a "
             "generic model failure.")


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def nc(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def orient_mult(lref, lalt, fref, falt):
    lref, lalt, fref, falt = (str(x).upper() for x in (lref, lalt, fref, falt))
    if (lref, lalt) == (fref, falt):
        return 1
    if (lref, lalt) == (falt, fref):
        return -1
    cf_ref, cf_alt = COMP.get(fref), COMP.get(falt)
    if (lref, lalt) == (cf_ref, cf_alt):
        return 1
    if (lref, lalt) == (cf_alt, cf_ref):
        return -1
    return 0


def auroc(y, s, ok=None):
    if ok is None:
        ok = ~np.isnan(s)
    if ok.sum() == 0 or len(np.unique(y[ok])) < 2:
        return np.nan
    return float(roc_auc_score(y[ok], s[ok]))


def block_bootstrap(y, score, groups, ok, nboot=NBOOT, seed=123):
    rng = np.random.default_rng(seed)
    yo, so, go = y[ok], score[ok], np.asarray(groups)[ok]
    keys = pd.unique(go)
    aucs = []
    for _ in range(nboot):
        pick = rng.choice(keys, size=len(keys), replace=True)
        idx = np.concatenate([np.where(go == c)[0] for c in pick])
        if len(np.unique(yo[idx])) < 2:
            continue
        aucs.append(roc_auc_score(yo[idx], so[idx]))
    if not aucs:
        return (np.nan, np.nan)
    return (float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)))


def perm_null_p(y, score, ok, obs, nperm=NPERM, seed=99):
    rng = np.random.default_rng(seed)
    yo, so = y[ok], score[ok]
    null = []
    for _ in range(nperm):
        yp = rng.permutation(yo)
        if len(np.unique(yp)) < 2:
            continue
        null.append(roc_auc_score(yp, so))
    null = np.asarray(null)
    a = max(obs, 1 - obs)
    hit = np.sum((null >= a) | (null <= 1 - a))
    return (1 + int(hit)) / (1 + len(null)), (float(np.mean(null)) if len(null) else np.nan), len(null)


def positional_baseline(df, y, groups):
    cols = [pd.to_numeric(df.get("maf"), errors="coerce"),
            pd.to_numeric(df.get("tss_or_peak_dist"), errors="coerce"),
            pd.to_numeric(df.get("tss_or_peak_dist"), errors="coerce").abs()]
    X = np.nan_to_num(np.column_stack([c.to_numpy() for c in cols]), nan=0.0)
    oof = np.full(len(df), np.nan)
    for held in pd.unique(groups):
        tr, te = groups != held, groups == held
        if tr.sum() < 5 or len(np.unique(y[tr])) < 2:
            continue
        m = Pipeline([("sc", StandardScaler()),
                      ("clf", LogisticRegression(max_iter=2000, class_weight="balanced"))])
        m.fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    return oof


def build():
    lab = pd.read_csv(LAB, sep="\t")
    lab["chr"] = lab["chr"].map(nc)
    tru = pd.read_csv(TRUTH, sep="\t")
    tru["chr"] = tru["chr"].map(nc)
    enr = lab.merge(
        tru[["variant_id", "signal_id", "locus", "study", "trait", "ancestry",
             "susie_cs", "strand_ambiguous"]].drop_duplicates("variant_id"),
        left_on="variant", right_on="variant_id", how="left")
    # re-assert: drop any strand-ambiguous (verified already 0)
    n0 = len(enr)
    enr = enr[enr["strand_ambiguous"] != True].reset_index(drop=True)  # noqa: E712
    log(f"strand-ambiguous dropped: {n0 - len(enr)} (expected 0)")

    bz = pd.read_csv(BZ, sep="\t")
    bz["chr"] = bz["chr"].map(nc)
    bz = bz.drop_duplicates("variant_id")
    m = enr.merge(bz, left_on="variant", right_on="variant_id",
                  how="left", suffixes=("", "_bz"))
    # orient each zero-shot signed delta into the label ref->alt frame
    om = np.array([orient_mult(r["ref"], r["alt"], r.get("hg38_ref"), r.get("hg38_alt"))
                   for _, r in m.iterrows()], dtype=float)
    m["_om"] = om
    for col in ZS_MODELS.values():
        m[col] = pd.to_numeric(m[col], errors="coerce") * m["_om"].replace(0, np.nan)
    # AlphaGenome native eQTL scorer — ZERO-SHOT annotation only (own orientation)
    ag = pd.read_csv(AG_FILE, sep="\t").drop_duplicates("variant_id")
    ag = ag.rename(columns={"hg38_ref": "ag_ref", "hg38_alt": "ag_alt"})
    m = m.merge(ag[["variant_id", "ag_ref", "ag_alt"] + list(AG_ZS.values())],
                left_on="variant", right_on="variant_id", how="left", suffixes=("", "_ag"))
    om_ag = np.array([orient_mult(r["ref"], r["alt"], r.get("ag_ref"), r.get("ag_alt"))
                      for _, r in m.iterrows()], dtype=float)
    m["_om_ag"] = om_ag
    for col in AG_ZS.values():
        m[col] = pd.to_numeric(m[col], errors="coerce") * pd.Series(om_ag, index=m.index).replace(0, np.nan)
    # canonical physical key (chr:pos:sorted-alleles) for the variant-level sensitivity
    def canon(r):
        a, b = sorted([str(r["ref"]).upper(), str(r["alt"]).upper()])
        return f"{r['chr']}:{int(r['pos_hg38'])}:{a}:{b}"
    m["phys_key"] = [canon(r) for _, r in m.iterrows()]
    return m


def run_level(df, level, cluster_col):
    y = (pd.to_numeric(df["sign"], errors="coerce").to_numpy() > 0).astype(int)
    chrom = df["chr"].astype(str).to_numpy()
    clust = df[cluster_col].astype(str).to_numpy()
    log(f"[{level}] n={len(df)} n_pos={int(y.sum())} n_chrom={df['chr'].nunique()} "
        f"n_cluster={pd.unique(clust).size}")
    rows = []
    scores = {}
    lic_map = {**{k: "open_zeroshot" for k in ZS_MODELS},
               **{k: "alphagenome_zeroshot_annotation" for k in AG_ZS}}
    for label, col in {**ZS_MODELS, **AG_ZS}.items():
        s = pd.to_numeric(df.get(col), errors="coerce").to_numpy()
        scores[label] = s
        ok = ~np.isnan(s)
        a = auroc(y, s, ok)
        lo, hi = block_bootstrap(y, s, chrom, ok)          # by chromosome
        clo, chi = block_bootstrap(y, s, clust, ok)        # by signal cluster
        p, nm, nn = perm_null_p(y, s, ok, a)
        rows.append(dict(model=label, auroc=a, ci_lo=lo, ci_hi=hi,
                         cluster_ci_lo=clo, cluster_ci_hi=chi,
                         perm_p=p, perm_null_mean=nm, n_perm=nn, n=int(ok.sum()),
                         meets_bar=bool(lo == lo and lo > PRIMARY_BAR), license=lic_map[label]))
    pos = positional_baseline(df, y, chrom)
    ok = ~np.isnan(pos)
    a = auroc(y, pos, ok)
    lo, hi = block_bootstrap(y, pos, chrom, ok)
    rows.append(dict(model="positional_MAF_tssdist_logistic", auroc=a, ci_lo=lo, ci_hi=hi,
                     cluster_ci_lo=np.nan, cluster_ci_hi=np.nan, perm_p=np.nan,
                     n=int(ok.sum()), meets_bar=False, license="baseline"))
    primary = [r for r in rows if r["model"] == "borzoi_logsed_liver"][0]
    return dict(level=level, n=len(df), n_pos=int(y.sum()),
                n_signal=int(pd.unique(clust).size), per_model=rows,
                primary_model="borzoi_logsed_liver", primary_auroc=primary["auroc"],
                primary_ci=[primary["ci_lo"], primary["ci_hi"]],
                primary_meets_bar=primary["meets_bar"])


def main():
    m = build()
    # signal level: one row per Broadaway signal_id (drop multi-gene/study replication)
    sig = m.sort_values("pvalue").drop_duplicates("signal_id").reset_index(drop=True)
    sig_level = run_level(sig, "signal_level_one_per_signal_id", "signal_id")
    # physical-variant sensitivity: one row per canonical physical SNP
    phys = m.sort_values("pvalue").drop_duplicates("phys_key").reset_index(drop=True)
    phys_level = run_level(phys, "physical_variant_sensitivity", "phys_key")
    # full (pseudo-replicated) for comparison to the old 211-row number
    full_level = run_level(m, "all_rows_pseudoreplicated_211", "signal_id")

    verdict = dict(
        script="85_eqtl_signal_diagnose.py",
        supersedes="eQTL arm of 80_direction_hardened.py (AG-trained, ToS-noncompliant)",
        tos_compliance=("Headline = zero-shot Borzoi logSED (open). AlphaGenome native eQTL scorer "
                        "included as zero-shot ANNOTATION only, never trained on. " + AG_NOTICE),
        framing=("Diagnostic, NOT a model verdict: in this Broadaway/colocalization-derived "
                 "liver benchmark, zero-shot Borzoi logSED sign-auROC did NOT meet the "
                 "prespecified 0.65 threshold. AlphaGenome's own native eQTL scorer, run zero-shot "
                 "on the SAME panel, is included to test whether the sub-threshold result reflects "
                 "panel ascertainment (independent liver-microarray marginal leads, out of the "
                 "fine-mapped-GTEx distribution where AG reports ~0.80) rather than model failure."),
        prespecified_bar=PRIMARY_BAR, nboot=NBOOT, nperm=NPERM,
        signal_level=sig_level, physical_variant=phys_level,
        all_rows_pseudoreplicated=full_level,
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(verdict, fh, indent=2)
    log(f"wrote {OUT_JSON}")
    flat = []
    for lv in (sig_level, phys_level, full_level):
        for r in lv["per_model"]:
            flat.append({"level": lv["level"], **r})
    pd.DataFrame(flat).to_csv(OUT_TSV, sep="\t", index=False)
    log(f"wrote {OUT_TSV}")
    for lv in (sig_level, phys_level, full_level):
        log(f"===== {lv['level']} (n={lv['n']} n_signal={lv['n_signal']} n_pos={lv['n_pos']}) =====")
        for r in lv["per_model"]:
            log(f"  {r['model']:34s} auROC={r['auroc']:.4f} "
                f"chrCI[{r['ci_lo']:.3f},{r['ci_hi']:.3f}] perm_p={r.get('perm_p', float('nan'))} n={r['n']}")
        log(f"  PRIMARY borzoi_logsed_liver={lv['primary_auroc']:.4f} "
            f"meets 0.65 bar={lv['primary_meets_bar']}")


if __name__ == "__main__":
    main()
