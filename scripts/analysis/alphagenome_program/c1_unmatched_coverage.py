#!/usr/bin/env python3
"""C1 supplementary: each model on its OWN maximal coverage, unmatched.

The C1 headline is the matched set. This script reports what each model scores when it is not
forced onto the intersection, so the cost of matching is visible and prediction C1.2 ("every
model's Spearman at the imbalance sites is below its caQTL value") can be evaluated per model
rather than only on the two-model matched set. Every row here is labelled unmatched and never
enters a contrast.

Reads only tables written by c1_matched_endogenous.py plus the endpoint-1 label assembly rebuilt
the same way, so the loader is identical to the one that produced the observed values.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import c1_matched_endogenous as C1


def per_model(df: pd.DataFrame, models: list[str], label_col: str, weight_col: str | None,
              tag: str, rows: list) -> None:
    for m in models:
        ok = df[df[m].notna() & df[label_col].notna()].copy()
        if len(ok) < 10:
            rows.append(dict(set=tag, model=m, n_variants=len(ok), n_blocks=np.nan,
                             signed_spearman=np.nan, sign_concordance=np.nan,
                             ci_lo=np.nan, ci_hi=np.nan, note="fewer than 10 scored variants"))
            continue
        ok = ok.sort_values("block_1mb", kind="mergesort").reset_index(drop=True)
        s = ok[m].values.astype(float)
        y = ok[label_col].values.astype(float)
        w = ok[weight_col].values.astype(float) if weight_col else None
        rho = C1.weighted_spearman(s, y, w) if w is not None else C1.signed_spearman(s, y)
        d = np.full(C1.N_BOOT, np.nan)
        bb = C1.BlockBootstrap(ok["block_1mb"].values, C1.N_BOOT, C1.BOOT_SEED)
        for i, idx in enumerate(bb.draws()):
            ww = w[idx] if w is not None else None
            d[i] = C1.weighted_spearman(s[idx], y[idx], ww) if ww is not None \
                else C1.signed_spearman(s[idx], y[idx])
        lo, hi = C1.ci(d)
        rows.append(dict(set=tag, model=m, n_variants=len(ok),
                         n_blocks=int(ok.block_1mb.nunique()),
                         signed_spearman=rho,
                         sign_concordance=C1.sign_concordance(s, y, w),
                         ci_lo=lo, ci_hi=hi, note=""))
        C1.log(f"[{tag}] {m}: n={len(ok)} blocks={ok.block_1mb.nunique()} "
               f"rho={rho:.4f} [{lo:.4f},{hi:.4f}]")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    tables = os.path.join(args.outdir, "tables")

    # endpoint 1: rebuild the full label assembly through the same loader
    cbp_a = C1.load_cbp(C1.CBP_ADULT, "chrombpnet_adult_hep")
    cbp_b = C1.load_cbp(C1.CBP_281367, "chrombpnet_adult_hep_gse281367")
    dirfeat = pd.read_csv(C1.DIRFEAT, sep="\t")
    dirfeat["key"] = dirfeat["canon"].astype(str)
    dirfeat = dirfeat.rename(columns={
        "cbp_logfc_labelframe": "chrombpnet_hepg2",
        "borzoi_atac_delta_labelframe": "borzoi_atac",
        "borzoi_dnase_delta_labelframe": "borzoi_dnase"})
    atlas = C1.load_atlas_chunks()
    lab1 = cbp_a[cbp_a.label == 1][["key", "chr_x", "pos_hg38", "ref", "alt", "beta_alt",
                                    "q_value", "chrombpnet_adult_hep"]].copy()
    lab1 = lab1.merge(cbp_b[["key", "chrombpnet_adult_hep_gse281367"]], on="key", how="left")
    lab1 = lab1.merge(dirfeat[["key", "chrombpnet_hepg2", "borzoi_atac", "borzoi_dnase"]],
                      on="key", how="left")
    lab1 = lab1.merge(atlas, on="key", how="left")
    lab1["block_1mb"] = [C1.block_of(c, p) for c, p in zip(lab1.chr_x, lab1.pos_hg38)]
    lab1["control_allele_identity"] = [
        C1.NT_RANK.get(str(a).upper(), np.nan) - C1.NT_RANK.get(str(r).upper(), np.nan)
        for r, a in zip(lab1.ref, lab1.alt)]
    lab1["control_position"] = (lab1.pos_hg38 % 1000) / 1000.0 - 0.5

    rows = []
    per_model(lab1, C1.MODELS_E1, "beta_alt", None, "endpoint1_currin_unmatched", rows)

    ase = pd.read_csv(os.path.join(tables, "endpoint2_ase_union_scores.tsv"), sep="\t")
    e2_models = [m for m in
                 ["chrombpnet_adult_hep", "chrombpnet_adult_hep_gse281367",
                  "alphagenome_atac_liver", "alphagenome_dnase_liver", "chrombpnet_hepg2",
                  "borzoi_atac", "borzoi_dnase", "control_allele_identity", "control_position"]
                 if m in ase.columns]
    per_model(ase, e2_models, "mean_log2_alt_over_ref", "precision_weight",
              "endpoint2_allelic_imbalance_unmatched", rows)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(tables, "c1_unmatched_per_model.tsv"), sep="\t", index=False)

    # q-value sanity: the endpoint-1 label set should be the FDR<5% lead set
    q = pd.to_numeric(lab1.q_value, errors="coerce")
    with open(os.path.join(args.outdir, "label_qvalue_check.txt"), "w") as fh:
        fh.write(f"n_leads\t{len(lab1)}\n")
        fh.write(f"q_max\t{float(q.max())}\n")
        fh.write(f"n_q_ge_0.05\t{int((q >= 0.05).sum())}\n")
        fh.write(f"n_q_missing\t{int(q.isna().sum())}\n")
    C1.log("done")


if __name__ == "__main__":
    main()
