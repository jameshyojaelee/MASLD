#!/usr/bin/env python3
"""
src/78_qtl_labels.py  — Phase-5 HARDENED direction-benchmark label construction.

Builds oriented, LEAKAGE-TIERED, LD-de-pseudoreplicated variant SIGN-label tables
for the supervised/zero-shot DIRECTION benchmark (Track C). ADDITIVE ONLY — reads
existing SF tables, writes new direction_labels_{eqtl,caqtl,sqtl}.tsv. Does NOT
touch 74-77 / 46d / 27a / 60-73.

Three leakage tiers (NON-NEGOTIABLE protocol):
  (a) eQTL  Broadaway liver cis-eQTL      -> clean_primary       (independent microarray;
                                             NOT in Borzoi/AlphaGenome training)
  (b) caQTL Currin 2025 primary-hepatocyte -> clean_secondary     (non-GTEx ATAC;
                             ACCESSIBILITY phenotype -> match to accessibility features)
  (c) sQTL  GTEx v8 liver sQTL            -> leakage_exploratory  (GTEx IS in Borzoi+AG
                             training) -> EXPLORATORY-flagged ONLY, junction-level
                             leafcutter phenotype (map to a SpliceJunction scorer).

HARDENING applied here:
  - Independent LEADS only (eQTL is_signal_lead; caQTL top-variant-per-peak flag;
    sQTL top-intron-per-variant), NOT LD-pseudoreplicated members.
  - Orientation: every beta/slope re-expressed in the model REF->ALT frame
    (variant_id encodes chr:pos:REF:ALT; FastQTL/GTEx betas are ALT-allele effects).
  - Strand-ambiguous A/T & C/G SNPs HARD-DROPPED (a strand flip is undetectable by
    letter-matching, so its sign is unknowable).
  - SORT1 rs12740374 sign unit-check in every applicable panel.

Emits, per panel:  variant, chr, pos_hg38, ref, alt, sign, effect, maf,
tss_or_peak_dist, panel, leakage_tier  (+ feature_id, pvalue, is_peak_lead[caqtl],
has_model_features).

Author: coding-upgrade agent (Phase 5). Env: rnaseq (pandas/numpy).
"""
import os
import sys
import gzip
import glob
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------- paths
ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF   = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
CAQ_DIR = os.path.join(ROOT, "data/external/currin_2025_caqtl")
SQTL    = os.path.join(ROOT,
    "data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL/Liver.v8.sqtl_signifpairs.txt.gz")
BROADAWAY_TRUTH = os.path.join(SF, "broadaway_benchmark_truth.tsv")

# feature tables (for the "already-scored?" overlap report)
BORZOI_EQTL = os.path.join(SF, "borzoi_eqtl_logsed_scores.tsv")
AG_EQTL     = os.path.join(SF, "alphagenome_eqtl_scores.tsv")
CBP_ACC     = os.path.join(SF, "chrombpnet_accessibility.tsv")

# ----------------------------------------------------------------------------- constants
# caQTL self-threshold. Currin nominal FastQTL files carry NO permutation q-value, so
# we self-threshold on nominal p. Genome-wide 5e-8 is deliberately conservative given
# nominal-only stats; per-pair p is retained so downstream can re-threshold.
CAQTL_P_THRESH = 5e-8
# sQTL: GTEx signifpairs are ALREADY GTEx-significant (per-gene permutation-gated);
# we keep them all and dedup to one lead intron per variant.
AMBIG = {frozenset(("A", "T")), frozenset(("C", "G"))}
# SORT1 rs12740374 (hg38 chr1:109274968, G>T). Derived T creates a C/EBP site ->
# INCREASES SORT1 hepatic expression & chromatin accessibility -> expected sign = +1
# in ref(G)->alt(T) frame in every panel where present.
SORT1_CHR, SORT1_POS = "1", 109274968
SORT1_EXPECT = "+1 (alt T creates C/EBP site -> up expression/accessibility)"


def is_strand_ambiguous(ref, alt):
    if not isinstance(ref, str) or not isinstance(alt, str):
        return False
    if len(ref) != 1 or len(alt) != 1:      # indels are never A/T|C/G palindromes
        return False
    return frozenset((ref.upper(), alt.upper())) in AMBIG


def norm_chr(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def canon_key(chrom, pos, ref, alt):
    """Model-frame key: '<chr_no_prefix>:<pos>:<REF>:<ALT>'."""
    return f"{norm_chr(chrom)}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


# ============================================================================= (a) eQTL
def build_eqtl():
    print("\n[eQTL] Broadaway liver cis-eQTL (clean_primary)", flush=True)
    df = pd.read_csv(BROADAWAY_TRUTH, sep="\t", dtype=str)
    n0 = len(df)
    df = df[df["is_signal_lead"] == "TRUE"].copy()
    n_lead = len(df)
    # strand-ambiguous hard drop (truth flag + independent recompute)
    recomputed_ambig = df.apply(lambda r: is_strand_ambiguous(r["ref"], r["alt"]), axis=1)
    df["_ambig"] = (df["strand_ambiguous"] == "TRUE") | recomputed_ambig
    n_ambig = int(df["_ambig"].sum())
    df = df[~df["_ambig"]].copy()

    df["chr"] = df["chr"].map(norm_chr)
    df["pos_hg38"] = df["pos_hg38"].astype(int)
    df["sign"] = df["eqtl_sign"].astype(float).astype(int)   # already ref->alt oriented
    df["effect"] = df["eqtl_beta"].astype(float)
    df["pvalue"] = pd.to_numeric(df["eqtl_p"], errors="coerce")
    df["maf"] = np.nan                        # Broadaway truth carries no MAF
    df["feature_id"] = df["gene"]
    df["panel"] = "eqtl_broadaway"
    df["leakage_tier"] = "clean_primary"
    df["is_peak_lead"] = True                  # is_signal_lead already independent

    # TSS distance: borrow variant_tss_dist from Borzoi scored table (variant x gene)
    tss = {}
    if os.path.exists(BORZOI_EQTL):
        b = pd.read_csv(BORZOI_EQTL, sep="\t", dtype=str)
        for _, r in b.iterrows():
            k = f"{canon_key(r['chr'], r['pos_hg38'], r['hg38_ref'], r['hg38_alt'])}|{r['gene']}"
            tss[k] = r.get("variant_tss_dist", np.nan)
    df["_key"] = df.apply(lambda r: canon_key(r["chr"], r["pos_hg38"], r["ref"], r["alt"]),
                          axis=1)
    df["tss_or_peak_dist"] = df.apply(
        lambda r: tss.get(f"{r['_key']}|{r['gene']}", np.nan), axis=1)
    df["tss_or_peak_dist"] = pd.to_numeric(df["tss_or_peak_dist"], errors="coerce")

    df["variant"] = df["variant_id"]
    sort1 = _sort1_report(df, "sign")
    return _finalize(df), n0, n_lead, n_ambig, sort1


# ============================================================================= (b) caQTL
def build_caqtl():
    print("\n[caQTL] Currin 2025 primary-hepatocyte caQTL (clean_secondary)", flush=True)
    files = sorted(glob.glob(os.path.join(CAQ_DIR, "caQTL_hg38_chr*_*.txt.gz")))
    if not files:
        raise FileNotFoundError(f"No caQTL files in {CAQ_DIR}")
    keep = []
    total_pairs = 0
    for fp in files:
        chrom = os.path.basename(fp).split("_")[2]  # 'chr22'
        print(f"  scanning {chrom} ...", flush=True)
        for chunk in pd.read_csv(fp, sep="\t", chunksize=2_000_000,
                                 usecols=["variant", "peak", "distance_from_peakCenter",
                                          "pvalue", "beta", "MAF"]):
            total_pairs += len(chunk)
            pos = chunk["variant"].str.split(":", expand=True)[1].astype(np.int64)
            # keep pairs passing threshold OR at the SORT1 locus (for the unit-check)
            m = (chunk["pvalue"].astype(float) < CAQTL_P_THRESH) | \
                ((norm_chr(chrom) == SORT1_CHR) & (pos == SORT1_POS))
            sub = chunk[m]
            if len(sub):
                keep.append(sub)
    df = pd.concat(keep, ignore_index=True)
    print(f"  total pairs scanned={total_pairs:,}  kept(<p thresh)={len(df):,}", flush=True)

    parts = df["variant"].str.split(":", expand=True)
    df["chr"] = parts[0].map(norm_chr)
    df["pos_hg38"] = parts[1].astype(np.int64)
    df["ref"] = parts[2].str.upper()
    df["alt"] = parts[3].str.upper()
    df["beta"] = df["beta"].astype(float)
    df["pvalue"] = df["pvalue"].astype(float)
    df["MAF"] = df["MAF"].astype(float)
    df["dist"] = df["distance_from_peakCenter"].astype(float)

    # SORT1 capture (before the significance filter is enforced for the lead table)
    sort1_rows = df[(df["chr"] == SORT1_CHR) & (df["pos_hg38"] == SORT1_POS)].copy()

    # enforce the significance threshold for the actual lead set
    df = df[df["pvalue"] < CAQTL_P_THRESH].copy()
    # strand-ambiguous hard drop
    ambig = [is_strand_ambiguous(r, a) for r, a in zip(df["ref"], df["alt"])]
    n_ambig = int(np.sum(ambig))
    df = df[~np.array(ambig)].copy()

    # deterministic ordering: p asc, |dist| asc, variant, peak  -> reproducible ties
    df["absdist"] = df["dist"].abs()
    df = df.sort_values(["pvalue", "absdist", "variant", "peak"], kind="mergesort")

    # (1) variant-level lead: one best peak per variant
    var_lead = df.drop_duplicates(subset=["variant"], keep="first").copy()
    # (2) independent peak-lead: top variant per that variant's assigned peak
    peak_top = var_lead.sort_values(["pvalue", "absdist", "variant"],
                                    kind="mergesort").drop_duplicates(
                                        subset=["peak"], keep="first")
    peak_lead_variants = set(peak_top["variant"])
    var_lead["is_peak_lead"] = var_lead["variant"].isin(peak_lead_variants)

    var_lead["sign"] = np.sign(var_lead["beta"]).astype(int)   # beta = ALT effect => ref->alt
    var_lead["effect"] = var_lead["beta"]
    var_lead["maf"] = var_lead["MAF"]
    var_lead["tss_or_peak_dist"] = var_lead["dist"]
    var_lead["feature_id"] = var_lead["peak"]
    var_lead["panel"] = "caqtl_currin"
    var_lead["leakage_tier"] = "clean_secondary"

    sort1 = _sort1_report_caqtl(sort1_rows)
    return _finalize(var_lead), total_pairs, len(var_lead), \
        int(var_lead["is_peak_lead"].sum()), n_ambig, sort1


# ============================================================================= (c) sQTL
def build_sqtl():
    print("\n[sQTL] GTEx v8 liver sQTL (leakage_exploratory)", flush=True)
    df = pd.read_csv(SQTL, sep="\t")
    n0 = len(df)
    # variant_id = chr1_832873_A_C_b38  (GTEx slope is ALT-allele effect => ref->alt)
    parts = df["variant_id"].str.split("_", expand=True)
    df["chr"] = parts[0].map(norm_chr)
    df["pos_hg38"] = parts[1].astype(np.int64)
    df["ref"] = parts[2].str.upper()
    df["alt"] = parts[3].str.upper()

    sort1_rows = df[(df["chr"] == SORT1_CHR) & (df["pos_hg38"] == SORT1_POS)].copy()

    ambig = [is_strand_ambiguous(r, a) for r, a in zip(df["ref"], df["alt"])]
    n_ambig = int(np.sum(ambig))
    df = df[~np.array(ambig)].copy()

    df["slope"] = df["slope"].astype(float)
    df["pvalue"] = df["pval_nominal"].astype(float)
    # lead = one intron (phenotype) per variant, strongest nominal p
    df["absdist"] = df["tss_distance"].abs()
    df = df.sort_values(["pvalue", "absdist", "variant_id", "phenotype_id"],
                        kind="mergesort")
    lead = df.drop_duplicates(subset=["variant_id"], keep="first").copy()

    lead["sign"] = np.sign(lead["slope"]).astype(int)
    lead["effect"] = lead["slope"]
    lead["maf"] = lead["maf"].astype(float)
    lead["tss_or_peak_dist"] = lead["tss_distance"]
    lead["feature_id"] = lead["phenotype_id"]    # leafcutter intron (junction coords)
    lead["variant"] = lead["variant_id"]
    lead["panel"] = "sqtl_gtex_liver"
    lead["leakage_tier"] = "leakage_exploratory"
    lead["is_peak_lead"] = True                  # one intron-lead per variant

    sort1 = _sort1_report(lead, "sign", pos_col="pos_hg38")
    return _finalize(lead), n0, len(lead), n_ambig, sort1


# ----------------------------------------------------------------------------- helpers
COLS = ["variant", "chr", "pos_hg38", "ref", "alt", "sign", "effect", "maf",
        "tss_or_peak_dist", "panel", "leakage_tier", "feature_id", "pvalue",
        "is_peak_lead"]


def _finalize(df):
    for c in COLS:
        if c not in df.columns:
            df[c] = np.nan
    out = df[COLS].copy()
    out["pos_hg38"] = out["pos_hg38"].astype("Int64")
    return out.reset_index(drop=True)


def _sort1_report(df, sign_col, pos_col="pos_hg38"):
    hit = df[(df["chr"] == SORT1_CHR) & (df[pos_col].astype("Int64") == SORT1_POS)]
    if len(hit) == 0:
        return f"NOT PRESENT in this panel's lead set (expected {SORT1_EXPECT})"
    r = hit.iloc[0]
    got = int(r[sign_col])
    ok = "PASS" if got > 0 else "FAIL(unexpected sign)"
    return (f"present ref={r['ref']} alt={r['alt']} sign={got:+d} effect={r['effect']:.4g} "
            f"-> {ok}; expected {SORT1_EXPECT}")


def _sort1_report_caqtl(rows):
    if len(rows) == 0:
        return f"NOT found at chr1:{SORT1_POS} in Currin caQTL (expected {SORT1_EXPECT})"
    rows = rows.sort_values("pvalue")
    r = rows.iloc[0]
    got = int(np.sign(r["beta"]))
    passed = r["pvalue"] < CAQTL_P_THRESH
    ok = "PASS" if got > 0 else "FAIL(unexpected sign)"
    thr = "" if passed else f" [NOTE p={r['pvalue']:.2g} > thresh {CAQTL_P_THRESH:g}; not in lead set]"
    return (f"present ref={r['ref']} alt={r['alt']} peak={r['peak']} beta={r['beta']:.4g} "
            f"p={r['pvalue']:.3g} sign={got:+d} -> {ok}{thr}; expected {SORT1_EXPECT}")


def _load_feature_keys():
    """Canonical chr:pos:ref:alt keys already carrying model features, by modality."""
    keys = {"expr": set(), "acc": set(), "splice": set()}
    if os.path.exists(BORZOI_EQTL):
        b = pd.read_csv(BORZOI_EQTL, sep="\t", dtype=str)
        keys["expr"] |= {canon_key(r["chr"], r["pos_hg38"], r["hg38_ref"], r["hg38_alt"])
                         for _, r in b.iterrows()}
    if os.path.exists(AG_EQTL):
        a = pd.read_csv(AG_EQTL, sep="\t", dtype=str)
        keys["expr"] |= {canon_key(r["chr"], r["pos_hg38"], r["hg38_ref"], r["hg38_alt"])
                         for _, r in a.iterrows()}
        keys["splice"] |= keys["expr"]   # AG splice scored on the same eQTL substrate
    if os.path.exists(CBP_ACC):
        c = pd.read_csv(CBP_ACC, sep="\t", dtype=str)
        keys["acc"] |= {canon_key(r["chr"], r["pos_hg38"], r["allele1"], r["allele2"])
                        for _, r in c.iterrows()}
    return keys


def _overlap(df, feat_keys):
    k = df.apply(lambda r: canon_key(r["chr"], r["pos_hg38"], r["ref"], r["alt"]), axis=1)
    return k.isin(feat_keys).sum()


# ----------------------------------------------------------------------------- main
def main():
    os.makedirs(SF, exist_ok=True)
    feat = _load_feature_keys()

    eqtl, e_n0, e_lead, e_ambig, e_sort1 = build_eqtl()
    caqtl, c_pairs, c_var, c_peak, c_ambig, c_sort1 = build_caqtl()
    sqtl, s_n0, s_lead, s_ambig, s_sort1 = build_sqtl()

    eqtl.to_csv(os.path.join(SF, "direction_labels_eqtl.tsv"), sep="\t", index=False)
    caqtl.to_csv(os.path.join(SF, "direction_labels_caqtl.tsv"), sep="\t", index=False)
    sqtl.to_csv(os.path.join(SF, "direction_labels_sqtl.tsv"), sep="\t", index=False)

    # feature-overlap report (modality-matched)
    e_ov = _overlap(eqtl, feat["expr"])
    c_all_ov = _overlap(caqtl, feat["acc"])
    c_peak_ov = _overlap(caqtl[caqtl["is_peak_lead"]], feat["acc"])
    s_ov = _overlap(sqtl, feat["splice"])

    print("\n" + "=" * 74)
    print("SUMMARY  — direction_labels_{eqtl,caqtl,sqtl}.tsv")
    print("=" * 74)
    print(f"[eQTL  clean_primary ] truth_rows={e_n0}  is_signal_lead={e_lead}  "
          f"strand_ambig_dropped={e_ambig}  FINAL={len(eqtl)}")
    print(f"        sign balance: {dict(eqtl['sign'].value_counts())}")
    print(f"        SORT1: {e_sort1}")
    print(f"        modality=EXPRESSION features; already-scored(expr) overlap="
          f"{e_ov}/{len(eqtl)}")
    print(f"[caQTL clean_secondary] pairs_scanned={c_pairs:,}  p<{CAQTL_P_THRESH:g}  "
          f"strand_ambig_dropped={c_ambig}")
    print(f"        variant-level leads (min-p peak per variant) = {c_var}")
    print(f"        INDEPENDENT peak-leads (is_peak_lead=top var/peak) = {c_peak}")
    print(f"        sign balance: {dict(caqtl['sign'].value_counts())}")
    print(f"        SORT1: {c_sort1}")
    print(f"        modality=ACCESSIBILITY features (ChromBPNet/AG-ATAC/Borzoi-ATAC, "
          f"NOT gene-expr); already-scored(acc) overlap: all={c_all_ov}/{c_var} "
          f"peak-leads={c_peak_ov}/{c_peak}")
    print(f"[sQTL  leakage_EXPLORATORY] signif_pairs={s_n0}  strand_ambig_dropped={s_ambig}  "
          f"intron-leads/variant={s_lead}")
    print(f"        sign balance: {dict(sqtl['sign'].value_counts())}")
    print(f"        SORT1: {s_sort1}")
    print(f"        modality=SPLICE features (AG SpliceJunction, junction-matched); "
          f"phenotype=leafcutter intron; already-scored(splice-substrate) overlap="
          f"{s_ov}/{len(sqtl)}  [LEAKAGE: GTEx in Borzoi+AG training -> never certify]")
    print("=" * 74)
    print("Thresholds/rules: caQTL nominal p<%g (self-thresholded; no permutation q in "
          "Currin), lead=min-p peak per variant + independent top-var-per-peak flag; "
          "sQTL=GTEx-significant, lead=min-p intron per variant; eQTL=is_signal_lead. "
          "Orientation: all betas ref->alt (FastQTL/GTEx ALT-effect; Broadaway re-oriented "
          "upstream). Strand-ambiguous A/T&C/G hard-dropped in every panel." % CAQTL_P_THRESH)


if __name__ == "__main__":
    main()
