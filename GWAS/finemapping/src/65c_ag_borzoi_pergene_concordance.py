#!/usr/bin/env python3
"""
65c_ag_borzoi_pergene_concordance.py  --  AlphaGenome x Borzoi per-gene two-model
concordance SIDECAR (Stage-2b, annotation/nomination ONLY).

One row per (variant, eGene) Broadaway lead (+ shared anchors) recording whether the
two INDEPENDENT sequence-to-function models call the same DIRECTION for the hg38
ref->alt substitution, each model's internal confidence, a TOP-CONFIDENCE tier flag,
and the measured eQTL sign where available:
  - AlphaGenome : ag_logsed = GeneMaskLFCScorer signed liver-RNA log2 LFC (ag_gene_lfc);
                  ag_quantile = gnomAD common-variant quantile (signed [-1,1]).
  - Borzoi      : borzoi_ensemble_logsed = the ON-RECIPE 4-fold (fwd+RC) GTEx-liver
                  tissue-matched logSED (`borzoi_logsed_liver`, job 18759567 --
                  supersedes single-fold); borzoi_fold_agreement = # of 4 folds
                  agreeing on sign (0-4; the ensemble has no gnomAD quantile).

Both models are scored strictly hg38 REF->ALT, so sign(ag_logsed) vs
sign(borzoi_ensemble_logsed) is a direct MODEL-vs-MODEL comparison (rows with matching
hg38 alleles). The measured Broadaway eqtl_sign is EFFECT-ALLELE-oriented (+1 => effect
allele raises expr); model_oriented_sign re-orients the (agreed) model direction to the
effect allele so measured_concordant is an apples-to-apples model-vs-truth check.

TOP-CONFIDENCE TIER (both_high_conf): borzoi 4/4-fold-unanimous AND |ag_quantile|>=0.9
AND sign_agree. Given the blanket eQTL-direction gate FAILed for BOTH models (AG 0.565 /
Borzoi-ensemble 0.559), this tier is the only subset where a directional nomination is
defensible -- both models are internally confident AND agree.

GUARDRAIL (team-lead sign-off): annotation/nomination-only. NOT wired into the
convergence atlas (46d/27a/78) or any scored channel (circularity with COLOC +
epigenomic). Additive only (script 65c); writes one sidecar TSV + a README.
"""
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SEQ = f"{ROOT}/GWAS/finemapping/results/seqfunc"
AG = f"{SEQ}/alphagenome_eqtl_scores.tsv"
BZ = f"{SEQ}/borzoi_eqtl_logsed_scores.tsv"          # ENSEMBLE (on-recipe) -- canonical
TRUTH = f"{SEQ}/broadaway_benchmark_truth.tsv"
OUT = f"{SEQ}/ag_borzoi_pergene_concordance.tsv"
OUT_README = f"{SEQ}/ag_borzoi_pergene_concordance.README.txt"

AG_Q_THRESH = 0.9      # |ag_quantile| threshold for high-confidence
BZ_FOLD_UNANIMOUS = 4  # 4/4 folds


def strip_ver(e):
    return str(e).split(".")[0] if e is not None else e


def sgn(x):
    if pd.isna(x):
        return np.nan
    return 1 if x > 0 else (-1 if x < 0 else 0)


ag = pd.read_csv(AG, sep="\t", dtype=str)
bz = pd.read_csv(BZ, sep="\t", dtype=str)
tr = pd.read_csv(TRUTH, sep="\t", dtype=str)
for df, cols in ((ag, ["ag_gene_lfc", "ag_gene_quantile"]),
                 (bz, ["borzoi_logsed_liver", "fold_sign_agreement"]),
                 (tr, ["pp4_best", "eqtl_sign"])):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

for df in (ag, bz, tr):
    df["key"] = df.variant_id + "|" + df.ensembl.map(strip_ver).fillna("")

lead_keys = set(tr[tr.is_signal_lead == "TRUE"]["key"])
anchor_keys = set(ag[ag.is_anchor == "TRUE"]["key"]) | set(bz[bz.is_anchor == "TRUE"]["key"])
keys = lead_keys | anchor_keys

tri = tr.drop_duplicates("key").set_index("key")
agi = ag.drop_duplicates("key").set_index("key")
bzi = bz.drop_duplicates("key").set_index("key")

rows = []
for key in sorted(keys):
    a = agi.loc[key] if key in agi.index else None
    b = bzi.loc[key] if key in bzi.index else None
    t = tri.loc[key] if key in tri.index else None
    src = a if a is not None else b
    if src is None:
        continue

    gene = src.get("gene", "")
    ensembl = src.get("ensembl", "")
    variant_id = src.get("variant_id", key.split("|")[0])
    is_anchor = str((a.get("is_anchor", "FALSE") if a is not None else b.get("is_anchor", "FALSE")))
    coloc_pp4 = (t["pp4_best"] if t is not None else np.nan)
    eqtl_sign = (t["eqtl_sign"] if t is not None else np.nan)
    effect_allele = (str(t["effect_allele"]) if t is not None and "effect_allele" in t else "")

    ag_lfc = a["ag_gene_lfc"] if a is not None else np.nan
    ag_q = a["ag_gene_quantile"] if a is not None else np.nan
    bz_lsd = b["borzoi_logsed_liver"] if b is not None else np.nan
    bz_fold = b["fold_sign_agreement"] if b is not None else np.nan
    hg38_ref = str((a.get("hg38_ref", "") if a is not None else b.get("hg38_ref", "")))
    hg38_alt = str((a.get("hg38_alt", "") if a is not None else b.get("hg38_alt", "")))

    ag_sign = sgn(ag_lfc)
    bz_sign = sgn(bz_lsd)

    # allele consistency (both scored hg38 ref->alt): require identical hg38 alleles
    note = ""
    alleles_ok = True
    if a is not None and b is not None:
        if (str(a.get("hg38_ref", "")) != str(b.get("hg38_ref", "")) or
                str(a.get("hg38_alt", "")) != str(b.get("hg38_alt", ""))):
            alleles_ok = False
            note = f"hg38_allele_mismatch(ag={a.get('hg38_ref')}/{a.get('hg38_alt')}," \
                   f"bz={b.get('hg38_ref')}/{b.get('hg38_alt')})"

    if not pd.isna(ag_sign) and not pd.isna(bz_sign) and alleles_ok:
        sign_agree = bool(ag_sign == bz_sign)
        sign_agree_str = "TRUE" if sign_agree else "FALSE"
    else:
        sign_agree = False
        sign_agree_str = "NA"
        if pd.isna(ag_sign):
            note = (note + ";" if note else "") + "ag_unscored"
        if pd.isna(bz_sign):
            note = (note + ";" if note else "") + "borzoi_unscored"

    # TOP-CONFIDENCE tier: both models internally confident AND agree on direction
    both_high_conf = bool(
        sign_agree and not pd.isna(ag_q) and not pd.isna(bz_fold)
        and abs(ag_q) >= AG_Q_THRESH and bz_fold >= BZ_FOLD_UNANIMOUS)

    # re-orient the (agreed) model direction to the effect allele, compare to measured
    model_oriented_sign = np.nan
    measured_concordant = ""
    if sign_agree and effect_allele and effect_allele in (hg38_alt, hg38_ref):
        oriented = ag_sign if effect_allele == hg38_alt else -ag_sign
        model_oriented_sign = int(oriented)
        if not pd.isna(eqtl_sign):
            measured_concordant = "TRUE" if oriented == sgn(eqtl_sign) else "FALSE"

    rows.append(dict(
        gene=gene, ensembl=ensembl, variant_id=variant_id,
        coloc_pp4=("" if pd.isna(coloc_pp4) else round(float(coloc_pp4), 4)),
        ag_logsed_sign=("" if pd.isna(ag_sign) else int(ag_sign)),
        borzoi_ens_logsed_sign=("" if pd.isna(bz_sign) else int(bz_sign)),
        sign_agree=sign_agree_str,
        ag_quantile=("" if pd.isna(ag_q) else round(float(ag_q), 6)),
        borzoi_fold_agreement=("" if pd.isna(bz_fold) else int(bz_fold)),
        both_high_conf=str(both_high_conf).upper(),
        eqtl_sign_if_measured=("" if pd.isna(eqtl_sign) else int(sgn(eqtl_sign))),
        effect_allele=effect_allele, hg38_ref=hg38_ref, hg38_alt=hg38_alt,
        ag_logsed=("" if pd.isna(ag_lfc) else round(float(ag_lfc), 6)),
        borzoi_ensemble_logsed=("" if pd.isna(bz_lsd) else round(float(bz_lsd), 6)),
        model_oriented_sign=("" if pd.isna(model_oriented_sign) else int(model_oriented_sign)),
        measured_concordant=measured_concordant,
        is_anchor=is_anchor, note=note,
    ))

cols = ["gene", "ensembl", "variant_id", "coloc_pp4",
        "ag_logsed_sign", "borzoi_ens_logsed_sign", "sign_agree",
        "ag_quantile", "borzoi_fold_agreement", "both_high_conf", "eqtl_sign_if_measured",
        "effect_allele", "hg38_ref", "hg38_alt", "ag_logsed", "borzoi_ensemble_logsed",
        "model_oriented_sign", "measured_concordant", "is_anchor", "note"]
out = pd.DataFrame(rows, columns=cols)
out.to_csv(OUT, sep="\t", index=False)

# ---------------------------------------------------------------- report
comp = out[out.sign_agree.isin(["TRUE", "FALSE"])]
n_comp = len(comp)
n_agree = int((comp.sign_agree == "TRUE").sum())
pct_agree = round(100 * n_agree / n_comp, 2) if n_comp else float("nan")

hc = out[out.both_high_conf == "TRUE"]
n_hc_rows = len(hc)
n_hc_genes = hc.gene.nunique()
hc_meas = hc[hc.measured_concordant.isin(["TRUE", "FALSE"])]
n_hc_meas = len(hc_meas)
n_hc_conc = int((hc_meas.measured_concordant == "TRUE").sum())
pct_hc_conc = round(100 * n_hc_conc / n_hc_meas, 2) if n_hc_meas else float("nan")

from scipy.stats import binomtest
hc_conc_p = (binomtest(n_hc_conc, n_hc_meas, 0.5, alternative="greater").pvalue
             if n_hc_meas else float("nan"))

print("=" * 66)
print(f"wrote {OUT}  ({len(out)} rows)")
print(f"comparable (both scored, alleles match)     : {n_comp}")
print(f"two-model SIGN AGREEMENT                     : {n_agree}/{n_comp} = {pct_agree}%")
print("-" * 66)
print(f"TOP-CONFIDENCE tier (both_high_conf==TRUE)  : {n_hc_rows} leads / {n_hc_genes} genes")
print(f"  ...with a measured Broadaway eQTL         : {n_hc_meas}")
print(f"  ...model direction MATCHES measured sign  : {n_hc_conc}/{n_hc_meas} = {pct_hc_conc}%  "
      f"(binom p vs 0.5 = {hc_conc_p:.3g})")
print("-" * 66)
anc = out[out.is_anchor == "TRUE"]
for _, r in anc.iterrows():
    print(f"  ANCHOR {r['gene']:9} ag_sign={r['ag_logsed_sign']} bz_sign={r['borzoi_ens_logsed_sign']} "
          f"sign_agree={r['sign_agree']} both_high_conf={r['both_high_conf']}")
print("=" * 66)

with open(OUT_README, "w") as f:
    f.write(
        "ag_borzoi_pergene_concordance.tsv -- AlphaGenome x Borzoi two-model per-gene\n"
        "direction concordance SIDECAR (ANNOTATION / NOMINATION ONLY).\n\n"
        "GUARDRAIL: NOT wired into the convergence atlas (46d/27a/78) or any scored\n"
        "channel. The seqfunc layer is annotation/nomination-only, never a scored\n"
        "convergence vote (circularity with COLOC + epigenomic). Atlas integration is a\n"
        "separate sign-off-gated stage.\n\n"
        "One row per (variant, eGene) Broadaway lead (+ shared anchors, is_anchor=TRUE).\n\n"
        "SPEC COLUMNS:\n"
        "  gene, ensembl, variant_id     target eGene + lead variant.\n"
        "  coloc_pp4                     Broadaway SuSiE/ABF best PP.H4 (annotation only).\n"
        "  ag_logsed_sign                sign of AlphaGenome GeneMaskLFCScorer liver log2\n"
        "                                LFC, hg38 REF->ALT (+1/-1/0).\n"
        "  borzoi_ens_logsed_sign        sign of Borzoi ON-RECIPE 4-fold GTEx-liver logSED\n"
        "                                (borzoi_logsed_liver), hg38 REF->ALT.\n"
        "  sign_agree                    TRUE/FALSE/NA -- both scored, hg38 alleles match,\n"
        "                                signs equal (MODEL-vs-MODEL, ref->alt space).\n"
        "  ag_quantile                   AG gnomAD common-variant quantile (signed [-1,1]).\n"
        "  borzoi_fold_agreement         # of 4 Borzoi folds agreeing on sign (0-4).\n"
        "  both_high_conf                TRUE iff borzoi_fold_agreement==4 AND\n"
        "                                |ag_quantile|>=0.9 AND sign_agree. THE TOP-CONFIDENCE\n"
        "                                TIER: the only subset where a directional nomination\n"
        "                                is defensible (the blanket eQTL-direction gate FAILed\n"
        "                                for BOTH models: AG auROC 0.565 / Borzoi-ens 0.559).\n"
        "  eqtl_sign_if_measured         measured Broadaway eQTL sign, EFFECT-ALLELE-oriented\n"
        "                                (+1 => effect allele raises expr); '' if no lead eQTL.\n\n"
        "AUDIT COLUMNS:\n"
        "  effect_allele, hg38_ref, hg38_alt, ag_logsed, borzoi_ensemble_logsed  provenance.\n"
        "  model_oriented_sign           the agreed model direction RE-ORIENTED to the effect\n"
        "                                allele (so it is comparable to eqtl_sign_if_measured).\n"
        "  measured_concordant           TRUE/FALSE -- model_oriented_sign == measured eqtl_sign\n"
        "                                (only where models agree AND an eQTL is measured).\n"
        "  is_anchor, note               provenance.\n\n"
        "NOTE ON ORIENTATION: ag_logsed_sign and borzoi_ens_logsed_sign are hg38 REF->ALT.\n"
        "eqtl_sign_if_measured is EFFECT-ALLELE-oriented. Do NOT compare them directly --\n"
        "use model_oriented_sign (already re-oriented) for the model-vs-measured check.\n"
    )
print(f"wrote {OUT_README}")
