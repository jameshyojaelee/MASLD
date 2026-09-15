#!/usr/bin/env python3
"""Step 09 (Package A): molecular-mechanism landscape per eligible signal.

No argmax across modalities. Every channel keeps its own coverage; flags record agreement
and ambiguity. Signed values only for scorers the Atlas declares signed.

Outputs (tables/):
  signal_channel_profiles.tsv   signal × channel family: coverage, signed/magnitude (raw and quantile), completeness
  signal_profiles.tsv           one row per signal with the ambiguity flags and consequence-class mass
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict

import pandas as pd

import lib_atlas as la

P = la.prespec()
TABLES = la.out_root() / "tables"
CONSEQ = la.PROJECT / "RNA-seq/results/coloc_variant_classes/cs_member_annotation.csv"
OPP_W = P["coverage"]["opposing_sign_min_weight"]

CHANNELS = {
    "expression": ["RNA_SEQ"],
    "expression_active": ["RNA_SEQ_ACTIVE"],
    "splice_site_usage": ["SPLICE_SITE_USAGE"],
    "splice_junctions": ["SPLICE_JUNCTIONS"],
    "splice_sites": ["SPLICE_SITES"],
    "polyadenylation": ["POLYADENYLATION"],
    "accessibility_atac": ["ATAC"],
    "accessibility_dnase": ["DNASE"],
    "histone": ["CHIP_HISTONE"],
    "tf_binding": ["CHIP_TF"],
    "cage": ["CAGE"],
    "procap": ["PROCAP"],
    "contact": ["CONTACT_MAPS"],
}
STRONG_Q, QUIET_Q = 0.9, 0.5


def main() -> None:
    signals = {s["signal_uid"]: s for s in la.read_tsv(TABLES / "eligible_signals.tsv")}
    meta = {r["scorer"]: r for r in la.read_tsv(TABLES / "scorer_metadata.tsv")}
    ps = pd.read_csv(TABLES / "posterior_summaries.tsv.gz", sep="\t", keep_default_na=False, na_values=["nan", "NaN", ""])
    avi = pd.read_csv(TABLES / "avi_records.tsv.gz", sep="\t").set_index("variant_uid")
    liver = pd.read_csv(TABLES / "liver_summaries.tsv.gz", sep="\t", keep_default_na=False)
    rna_liver = liver[liver["scorer"] == "RNA_SEQ"]

    # weights (mapped, queried) per signal and consequence classes by hg19 key
    w_q = defaultdict(dict)
    src_of = {}
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["mapping_status"] == "mapped" and r["in_query_set"] == "True":
                w_q[r["signal_uid"]][r["variant_uid"]] = float(r["weight"])
            src_of[(r["signal_uid"], r["variant_uid"])] = r["source_variant_id"]
    conseq = {}
    with open(CONSEQ, newline="") as handle:
        for r in csv.DictReader(handle):
            conseq[(r["variant_key"], r["allele1"].upper(), r["allele2"].upper())] = r["coarse_class"]
            conseq[(r["variant_key"], r["allele2"].upper(), r["allele1"].upper())] = r["coarse_class"]

    channel_rows, profile_rows = [], []
    for sig, s in signals.items():
        sub = ps[ps["signal_uid"] == sig]
        target = s["ensembl"]
        prof = {"signal_uid": sig, "universe": s["universe"], "gwas_name": s["gwas_name"], "trait": s["trait"], "trait_class": s["trait_class"],
                "gene": s["gene"], "ensembl": target, "analysis_block": s["analysis_block"], "posterior_definition": s["posterior_definition"]}
        chan_q = {}
        for chan, scorers in CHANNELS.items():
            for scorer in scorers:
                if scorer not in meta:
                    continue
                sel = sub[(sub["scorer"] == scorer) & (sub["gene_scope"].isin(["target_gene", "region"]))]
                m = meta[scorer]
                liver_avail = int(m["n_primary_liver"]) + int(m["n_hepatocyte"]) > 0
                row = {"signal_uid": sig, "universe": s["universe"], "channel": chan, "scorer": scorer, "is_signed": m["is_signed"] == "True",
                       "liver_tracks_available": liver_avail, "surrogate_only": (not liver_avail) and int(m["n_hepg2"]) > 0,
                       "gene_scope": "", "C_liver": math.nan, "S_raw": math.nan, "A_raw": math.nan, "S_over_C": math.nan, "A_over_C": math.nan,
                       "S_quantile": math.nan, "A_quantile": math.nan, "A_quantile_over_C": math.nan, "complete_0.95": False,
                       "top_variant": "", "top_weight": math.nan, "state": "no_liver_track" if not liver_avail else "no_returned_row"}
                if len(sel):
                    r = sel.iloc[0]
                    c = float(r["C_liver"]) if r["C_liver"] != "" else 0.0
                    row.update(gene_scope=r["gene_scope"], C_liver=c, S_raw=_f(r["S_liver_median"]), A_raw=_f(r["A_liver_median"]),
                               S_over_C=_f(r["S_over_C"]), A_over_C=_f(r["A_over_C"]), S_quantile=_f(r["S_quantile_weighted"]),
                               A_quantile=_f(r["A_quantile_weighted"]), A_quantile_over_C=(_f(r["A_quantile_weighted"]) / c) if c > 0 else math.nan,
                               **{"complete_0.95": str(r["complete_0.95"]) == "True"}, top_variant=r["top_variant"], top_weight=_f(r["top_weight"]),
                               state="scored" if c > 0 else "returned_no_liver_score")
                channel_rows.append(row)
                chan_q[chan] = row
        # AVI, posterior weighted over queried variants that returned AVI
        wq = w_q[sig]
        avi_vals = {v: avi.at[v, "AVI_SCORE"] for v in wq if v in avi.index}
        avi_q = {v: avi.at[v, "AVI_SCORE_quantile"] for v in wq if v in avi.index}
        a_sum = la.posterior_summary(wq, avi_vals, is_signed=False)
        aq_sum = la.posterior_summary(wq, avi_q, is_signed=False)
        top = max(avi_vals, key=lambda v: abs(avi_vals[v])) if avi_vals else ""
        top_fi = ""
        if top:
            fi = avi.loc[top]
            feats = [c for c in avi.columns if c not in ("AVI_SCORE", "AVI_SCORE_quantile") and not c.startswith("MODEL_")]
            fi_cols = [c for c in avi.columns if c in ("MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS", "MAX_ABS_DNASE", "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CAGE", "MAX_ABS_PROCAP", "MAX_ABS_RNA_SEQ", "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE", "CACTUS_241_WAY", "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST", "PHASTCONS_470_WAY", "IS_INSERTION", "IS_DELETION")]
            ranked = sorted(((abs(float(fi[c])), c) for c in fi_cols if not math.isnan(float(fi[c]))), reverse=True)[:3]
            top_fi = ";".join(f"{c}={v:.3g}" for v, c in ranked)
        # consequence-class mass (credible-set annotation covers only part of the posterior)
        coding = noncoding = c_conseq = 0.0
        for v, w in wq.items():
            src = src_of.get((sig, v), "")
            parts = src.split(":")  # hg19:chr:pos:a1:a2 or hg19:chrN:pos:a1:a2
            if len(parts) >= 5:
                key = (f"{parts[1].replace('chr', '')}:{parts[2]}", parts[3].upper(), parts[4].upper())
                cls = conseq.get(key)
                if cls is not None:
                    c_conseq += w
                    if cls == "coding":
                        coding += w
                    else:
                        noncoding += w
        # opposing sign among weighty variants on the target gene (expression channel)
        opposing, n_pos, n_neg = False, 0, 0
        if target:
            g = rna_liver[(rna_liver["gene_id"] == target) & (rna_liver["variant_uid"].isin([v for v, w in wq.items() if w >= OPP_W]))]
            strong = g[g["liver_median_quantile"].abs() >= QUIET_Q]
            n_pos, n_neg = int((strong["liver_median_raw"] > 0).sum()), int((strong["liver_median_raw"] < 0).sum())
            opposing = n_pos > 0 and n_neg > 0
        # multiple candidate genes: genes with posterior-weighted |quantile|/C ≥ STRONG_Q on expression
        genes_strong = []
        for _, r in sub[(sub["scorer"] == "RNA_SEQ") & (sub["gene_scope"].isin(["target_gene", "other_gene"]))].iterrows():
            c = _f(r["C_liver"]); aq = _f(r["A_quantile_weighted"])
            if c > 0 and aq / c >= STRONG_Q:
                genes_strong.append(r["gene_id"])
        expr = chan_q.get("expression", {})
        splice_max = max([chan_q.get(k, {}).get("A_quantile_over_C", math.nan) for k in ("splice_site_usage", "splice_junctions", "polyadenylation")] + [math.nan], key=lambda x: -1 if math.isnan(x) else x)
        expr_q = expr.get("A_quantile_over_C", math.nan)
        liver_max = max([r.get("A_quantile_over_C", math.nan) for r in chan_q.values() if r.get("liver_tracks_available")] + [math.nan], key=lambda x: -1 if math.isnan(x) else x)
        prof.update({
            "n_queried_variants": len(wq), "queried_mass": sum(wq.values()),
            "avi_C": a_sum["coverage"], "avi_magnitude_weighted": a_sum["magnitude"], "avi_quantile_weighted": aq_sum["magnitude"],
            "avi_top_variant": top, "avi_top_score": avi_vals.get(top, math.nan), "avi_top_quantile": avi_q.get(top, math.nan), "avi_top_features": top_fi,
            "consequence_C": c_conseq, "coding_mass": coding, "noncoding_mass": noncoding,
            "expression_target_C": expr.get("C_liver", math.nan), "expression_target_S_over_C": expr.get("S_over_C", math.nan),
            "expression_target_A_quantile_over_C": expr_q, "expression_complete_0.95": expr.get("complete_0.95", False),
            "processing_max_A_quantile_over_C": splice_max,
            "flag_multiple_candidate_genes": len(genes_strong) > 1 or (len(genes_strong) == 1 and target and genes_strong[0] != target),
            "strong_expression_genes": ";".join(genes_strong),
            "flag_opposing_signs_target_expression": opposing, "n_positive_strong": n_pos, "n_negative_strong": n_neg,
            "flag_quiet_expression_with_processing": (not math.isnan(expr_q) and expr_q < QUIET_Q) and (not math.isnan(splice_max) and splice_max >= STRONG_Q),
            "flag_strong_avi_little_liver": (avi_q.get(top, 0) >= STRONG_Q) and (math.isnan(liver_max) or liver_max < QUIET_Q),
            "flag_surrogate_only_channels": ";".join(k for k, r in chan_q.items() if r.get("surrogate_only")),
        })
        profile_rows.append(prof)
    la.write_tsv_once(TABLES / "signal_channel_profiles.tsv", channel_rows, list(channel_rows[0].keys()))
    la.write_tsv_once(TABLES / "signal_profiles.tsv", profile_rows, list(profile_rows[0].keys()))
    la.log(f"Package A: {len(profile_rows)} signals, {len(channel_rows)} channel rows")


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


if __name__ == "__main__":
    main()
