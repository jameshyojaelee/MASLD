#!/usr/bin/env python3
"""
68_chrombpnet_accessibility.py  --  Stage-3 (Axis-3) seqfunc assembler.

Takes the kundajelab/variant-scorer per-fold outputs (already summarized across the
5 ENCODE HepG2 ATAC ChromBPNet folds and annotated with peak overlap + TF-motif-hit
overlap) and emits the canonical Axis-3 deliverable:

    chrombpnet_accessibility.tsv
    cols: variant, gene, cbp_abs_logfc, cbp_active_quantile, cbp_jsd,
          disrupted_tf_motif, in_hepg2_peak
    (+ carried: cbp_logfc[signed], cbp_abs_logfc_x_jsd[IES], abs_logfc.pval,
       jsd.pval, max_pip, var_class, variant_id_hg19, chr, pos_hg38,
       allele1, allele2, n_folds)

FRAMING (do NOT drop): the eQTL-DIRECTION gate FAILED for both Borzoi and
AlphaGenome (auROC ~0.53-0.56) and Borzoi ATAC accessibility direction also failed
(auROC 0.537). This axis is therefore STRICTLY direction-independent: it nominates
MAGNITUDE ("this variant disrupts a hepatocyte-accessible element") + WHICH-TF-
FOOTPRINT mechanism (disrupted_tf_motif) + WHICH-GENE (nearest effector carried from
the fine-mapping substrate). The signed cbp_logfc is carried for completeness ONLY;
it is NOT a validated direction call and must not be wired into a scored/convergence
channel.

ADDITIVE ONLY: reads the variant-scorer outputs + the fine-mapping substrate genemap;
writes one new TSV. Does not touch 46d/78/27a or 60-65.
"""
import argparse
import os
import sys
import pandas as pd
import numpy as np


def pick(df, *cands):
    for c in cands:
        if c in df.columns:
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations", required=True,
                    help="variant_annotation.py .annotations.tsv (mean-across-folds scores + peak/hit overlap)")
    ap.add_argument("--genemap", required=True,
                    help="variant_genemap.tsv sidecar (variant_id -> gene, max_pip, var_class, variant_id_hg19)")
    ap.add_argument("--n_folds", type=int, default=5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ann = pd.read_table(args.annotations)
    gm = pd.read_table(args.genemap)
    print("annotations shape:", ann.shape, "| cols:", list(ann.columns))

    # ---- map summarized-score columns (tolerant to naming) --------------------
    c_abs   = pick(ann, "abs_logfc.mean", "abs_logfc")
    c_logfc = pick(ann, "logfc.mean", "logfc")
    c_quant = pick(ann, "active_allele_quantile.mean", "active_allele_quantile")
    c_jsd   = pick(ann, "jsd.mean", "jsd")
    c_ies   = pick(ann, "abs_logfc_x_jsd.mean", "abs_logfc_x_jsd")
    c_abspv = pick(ann, "abs_logfc.mean.pval", "abs_logfc.pval")
    c_jsdpv = pick(ann, "jsd.mean.pval", "jsd.pval")
    for nm, c in [("abs_logfc", c_abs), ("active_allele_quantile", c_quant), ("jsd", c_jsd)]:
        if c is None:
            sys.exit(f"FATAL: required score column for '{nm}' not found in annotations.")

    # ---- peak overlap + motif hit --------------------------------------------
    c_peak = pick(ann, "peak_overlap")
    c_hitm = pick(ann, "hits_motifs")

    out = pd.DataFrame()
    out["variant"] = ann["variant_id"]
    out["chr"] = ann["chr"]
    out["pos_hg38"] = ann["pos"]
    out["allele1"] = ann["allele1"]
    out["allele2"] = ann["allele2"]
    out["cbp_abs_logfc"] = ann[c_abs]
    out["cbp_active_quantile"] = ann[c_quant]
    out["cbp_jsd"] = ann[c_jsd]
    out["cbp_logfc"] = ann[c_logfc] if c_logfc else np.nan            # signed; direction-INDEPENDENT, carried only
    out["cbp_abs_logfc_x_jsd"] = ann[c_ies] if c_ies else out["cbp_abs_logfc"] * out["cbp_jsd"]
    out["cbp_abs_logfc_pval"] = ann[c_abspv] if c_abspv else np.nan
    out["cbp_jsd_pval"] = ann[c_jsdpv] if c_jsdpv else np.nan
    out["in_hepg2_peak"] = ann[c_peak].astype(bool) if c_peak else False
    if c_hitm:
        dm = ann[c_hitm].fillna("-").astype(str)
        out["disrupted_tf_motif"] = dm.where(dm != "", "-")
    else:
        out["disrupted_tf_motif"] = "-"
    out["n_folds"] = args.n_folds

    # ---- merge gene / provenance from substrate genemap -----------------------
    gm2 = gm.rename(columns={"variant_id": "variant"})[
        ["variant", "gene", "max_pip", "var_class", "variant_id_hg19"]]
    out = out.merge(gm2, on="variant", how="left")

    # canonical column order (deliverable cols first)
    lead = ["variant", "gene", "cbp_abs_logfc", "cbp_active_quantile", "cbp_jsd",
            "disrupted_tf_motif", "in_hepg2_peak"]
    rest = [c for c in out.columns if c not in lead]
    out = out[lead + rest]
    out.sort_values("cbp_abs_logfc", ascending=False, inplace=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    # ---- report ---------------------------------------------------------------
    n = len(out)
    n_peak = int(out["in_hepg2_peak"].sum())
    n_motif = int((out["disrupted_tf_motif"] != "-").sum())
    print("=" * 64)
    print(f"WROTE {args.out}")
    print(f"n variants scored: {n}")
    print(f"n in a HepG2 ATAC peak (in_hepg2_peak): {n_peak} ({100*n_peak/max(n,1):.1f}%)")
    print(f"n overlapping a TF-motif instance (disrupted_tf_motif != '-'): {n_motif}")
    print(f"cbp_abs_logfc  median={out['cbp_abs_logfc'].median():.4f}  max={out['cbp_abs_logfc'].max():.4f}")
    print(f"cbp_jsd        median={out['cbp_jsd'].median():.4f}  max={out['cbp_jsd'].max():.4f}")
    # top footprint disruptions among in-peak variants (the meaningful set)
    inpk = out[out["in_hepg2_peak"] & (out["disrupted_tf_motif"] != "-")].copy()
    inpk = inpk.sort_values("cbp_abs_logfc_x_jsd", ascending=False)
    print("\nTop 20 TF-footprint disruptions (in-peak, ranked by |logFC|*JSD):")
    cols = ["variant", "gene", "disrupted_tf_motif", "cbp_abs_logfc", "cbp_jsd",
            "cbp_active_quantile", "cbp_abs_logfc_x_jsd", "max_pip"]
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(inpk[cols].head(20).to_string(index=False))
    # motif-frequency table among disrupted in-peak variants
    from collections import Counter
    mc = Counter()
    for s in inpk["disrupted_tf_motif"]:
        for m in str(s).split(","):
            if m and m != "-":
                mc[m] += 1
    print("\nMost frequently disrupted TF motifs (in-peak variants):")
    for m, c in mc.most_common(20):
        print(f"  {m}: {c}")
    print("=" * 64)


if __name__ == "__main__":
    main()
