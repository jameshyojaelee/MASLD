#!/usr/bin/env python
"""
67_pops_whichgene.py  --  PoPS Axis-2 orthogonal which-gene nomination
======================================================================
ADDITIVE (new src/67). NOMINATION / mechanism-class only -- NOT wired into the
convergence atlas or any scored channel (PoPS is orthogonal by construction:
genome-wide polygenic + gene features, no local eQTL / ATAC / sequence, so it
does not double-count the COLOC / epigenomic / sequence channels).

For every eQTL-absent fine-mapped locus (credible-set signal whose effector gene
does NOT colocalize; substrate colocalizes==FALSE regulatory variants), report:
  * pops_top_gene : PoPS-prioritised gene within +/-500 kb of the credible variant
  * nearest_gene  : nearest TSS in the PoPS gene universe (Weeks 2023 baseline)
  * agree         : PoPS-top == nearest  (the Weeks-recommended high-confidence combo)

Inputs
  --preds        MVP_NAFLD_EUR.preds        (pops.py output: ENSGID, PoPS_Score)
  --gene_annot   gene_annot_jun10.txt       (ENSGID NAME CHR START END TSS, GRCh37)
  --substrate    variant_substrate_hg38.tsv (carries pos_hg19, max_pip, var_class,
                                             colocalizes, assigned gene)
Output
  results/seqfunc/pops_nominations.tsv
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd

CLUMP_KB  = 500   # independent-locus clumping radius (kb) around the top-PIP seed
WINDOW_KB = 500   # PoPS candidate-gene window (kb) around the credible variant
HEADLINE_PIP = 0.10   # "genuinely fine-mapped" tier for headline stats


def log(m): print(f"[67_pops_whichgene] {m}", flush=True)


def clump(df, kb):
    """Greedy PIP-seeded clumping into independent loci (per chromosome).
    df: columns chr(int), pos(int), max_pip(float), + passthrough. Returns list
    of dicts: lead row + member positions. Highest-PIP variant seeds each locus."""
    loci = []
    for c, g in df.groupby("chr"):
        g = g.sort_values("max_pip", ascending=False).reset_index(drop=True)
        lead_pos = []
        assigned = [-1] * len(g)
        for i, row in g.iterrows():
            p = row["pos"]
            hit = -1
            for li, lp in enumerate(lead_pos):
                if abs(p - lp) <= kb * 1000:
                    hit = li; break
            if hit == -1:
                lead_pos.append(p)
                assigned[i] = len(lead_pos) - 1
            else:
                assigned[i] = hit
        g["_locus"] = assigned
        for li, gg in g.groupby("_locus"):
            seed = gg.iloc[0]  # highest PIP (sorted desc)
            loci.append({
                "chr": int(c),
                "pos_hg19": int(seed["pos"]),
                "credible_variant": seed["variant_id_hg19"],
                "lead_pip": float(seed["max_pip"]),
                "n_members": int(len(gg)),
                "substrate_assigned_gene": seed.get("gene", ""),
            })
    return loci


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preds", required=True)
    ap.add_argument("--gene_annot", required=True)
    ap.add_argument("--substrate", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # --- PoPS scores ---------------------------------------------------------
    preds = pd.read_csv(a.preds, sep="\t")
    preds = preds[["ENSGID", "PoPS_Score"]].dropna()
    preds["pctl"] = preds["PoPS_Score"].rank(pct=True)
    score = dict(zip(preds["ENSGID"], preds["PoPS_Score"]))
    pctl  = dict(zip(preds["ENSGID"], preds["pctl"]))
    log(f"PoPS scores: {len(score):,} genes")

    # --- gene annotation (GRCh37, PoPS universe) -----------------------------
    ga = pd.read_csv(a.gene_annot, sep="\t", dtype={"CHR": str})
    ga = ga[ga["ENSGID"].isin(score)].copy()
    ga["CHR"] = pd.to_numeric(ga["CHR"], errors="coerce")
    ga = ga.dropna(subset=["CHR", "TSS"])
    ga["CHR"] = ga["CHR"].astype(int)
    ga["TSS"] = ga["TSS"].astype(int)
    ga["PoPS_Score"] = ga["ENSGID"].map(score)
    # per-chr arrays for fast nearest / window queries
    by_chr = {c: g.sort_values("TSS").reset_index(drop=True) for c, g in ga.groupby("CHR")}
    log(f"PoPS gene universe with coords: {len(ga):,} genes (autosomes only)")

    # --- substrate: eQTL-absent regulatory -----------------------------------
    sub = pd.read_csv(a.substrate, sep="\t", dtype=str)
    sub["colocalizes_b"] = sub["colocalizes"].str.upper() == "TRUE"
    ea = sub[(~sub["colocalizes_b"]) & (sub["var_class"] == "regulatory")].copy()
    ea["chr"] = pd.to_numeric(ea["chr"].str.replace("chr", "", regex=False), errors="coerce")
    ea = ea[ea["chr"].notna()].copy()
    ea["chr"] = ea["chr"].astype(int)
    ea = ea[ea["chr"].between(1, 22)]                     # autosomes (PoPS universe)
    ea["pos"] = pd.to_numeric(ea["pos_hg19"], errors="coerce")
    ea["max_pip"] = pd.to_numeric(ea["max_pip"], errors="coerce").fillna(0.0)
    ea = ea.dropna(subset=["pos"]).copy()
    ea["pos"] = ea["pos"].astype(int)
    n_chrx = ((pd.to_numeric(sub["chr"].str.replace("chr", "", regex=False),
                             errors="coerce").isna()) &
              (~(sub["colocalizes"].str.upper() == "TRUE")) &
              (sub["var_class"] == "regulatory")).sum()
    log(f"eQTL-absent regulatory variants (autosomal): {len(ea):,}   "
        f"(non-autosomal excluded: {int(n_chrx)})")

    # --- clump into independent loci -----------------------------------------
    loci = clump(ea[["chr", "pos", "max_pip", "variant_id_hg19", "gene"]], CLUMP_KB)
    log(f"independent eQTL-absent loci (clumped @ {CLUMP_KB}kb): {len(loci):,}")

    # --- nominate ------------------------------------------------------------
    rows = []
    for L in loci:
        c, p = L["chr"], L["pos_hg19"]
        g = by_chr.get(c)
        if g is None or len(g) == 0:
            continue
        d = (g["TSS"] - p).abs()
        # nearest TSS in PoPS universe
        ni = d.idxmin()
        nearest_ens, nearest_name = g.loc[ni, "ENSGID"], g.loc[ni, "NAME"]
        nearest_dist = int(d.loc[ni])
        # candidate genes within window
        win = g[d <= WINDOW_KB * 1000]
        window_empty = len(win) == 0
        if window_empty:
            cand = g.loc[[ni]]            # fall back to nearest
        else:
            cand = win
        ti = cand["PoPS_Score"].idxmax()
        top_ens, top_name = cand.loc[ti, "ENSGID"], cand.loc[ti, "NAME"]
        top_score = float(cand.loc[ti, "PoPS_Score"])
        top_dist = int(abs(int(cand.loc[ti, "TSS"]) - p))
        agree = bool(top_ens == nearest_ens)
        rows.append({
            "locus": f"chr{c}:{p}",
            "credible_variant": L["credible_variant"],
            "pops_top_gene": top_name,
            "pops_score": round(top_score, 6),
            "nearest_gene": nearest_name,
            "agree": agree,
            "pops_and_nearest": top_name if agree else f"{top_name}|{nearest_name}",
            "chr": c,
            "pos_hg19": p,
            "lead_pip": round(L["lead_pip"], 4),
            "pops_top_ensgid": top_ens,
            "pops_top_pctl": round(float(pctl.get(top_ens, np.nan)), 4),
            "pops_top_dist_kb": round(top_dist / 1000, 1),
            "nearest_ensgid": nearest_ens,
            "nearest_dist_kb": round(nearest_dist / 1000, 1),
            "n_genes_in_window": int(len(win)),
            "n_locus_members": L["n_members"],
            "substrate_assigned_gene": L["substrate_assigned_gene"],
            "window_empty": window_empty,
        })

    out = pd.DataFrame(rows).sort_values(["lead_pip", "pops_score"],
                                         ascending=False).reset_index(drop=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)
    log(f"wrote {a.out}  ({len(out):,} loci)")

    # --- summary -------------------------------------------------------------
    def agr(df):
        return f"{df['agree'].mean()*100:.1f}% ({int(df['agree'].sum())}/{len(df)})" if len(df) else "n/a"
    hi = out[out["lead_pip"] >= HEADLINE_PIP]
    hi5 = out[out["lead_pip"] >= 0.5]
    log("================ SUMMARY ================")
    log(f"eQTL-absent loci with a PoPS gene       : {len(out):,}")
    log(f"  with non-empty +/-{WINDOW_KB}kb window     : {int((~out['window_empty']).sum()):,}")
    log(f"PoPS-vs-nearest agreement (all loci)     : {agr(out)}")
    log(f"  headline tier lead_pip>={HEADLINE_PIP} (n={len(hi)}) : {agr(hi)}")
    log(f"  strict   tier lead_pip>=0.5  (n={len(hi5)}) : {agr(hi5)}")
    log("========================================")


if __name__ == "__main__":
    main()
