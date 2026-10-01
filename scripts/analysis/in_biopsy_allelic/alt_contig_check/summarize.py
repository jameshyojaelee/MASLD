#!/usr/bin/env python3
"""Summarize the ALT/patch-contig undercount check.

For primary-chromosome genes, compares featureCounts with multimappers
(-M --fraction) to the default (multimappers dropped). Genes whose gene_name
also occurs on a non-primary contig ("ALT-copy") are compared with all other
primary genes, and both are joined to the canonical DEG table
(padj < 0.05 and |logFC| > 0.5, the definition in load_figure_data.R).
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

PRIMARY = {f"chr{i}" for i in list(range(1, 23)) + ["X", "Y", "M"]}


def genes(gtf):
    rows = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            a = f[8]
            rows.append((a.split('gene_id "', 1)[1].split('"', 1)[0], a.split('gene_name "', 1)[1].split('"', 1)[0],
                         a.split('gene_type "', 1)[1].split('"', 1)[0], f[0]))
    return pd.DataFrame(rows, columns=["gene_id", "gene_name", "gene_type", "seqname"])


def read_fc(path):
    d = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    return d.iloc[:, 5:].sum(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--counts-dir", required=True)
    ap.add_argument("--degs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-count", type=float, default=50)
    a = ap.parse_args()
    out = Path(a.out)

    g = genes(a.gtf)
    alt_names = set(g.loc[~g["seqname"].isin(PRIMARY), "gene_name"])
    prim = g[g["seqname"].isin(PRIMARY)].set_index("gene_id")
    prim["alt_copy"] = prim["gene_name"].isin(alt_names)

    per = []
    for f in sorted(Path(a.counts_dir).glob("*.default.txt")):
        c = f.name.split(".")[0]
        d, m = read_fc(f), read_fc(f.with_name(f"{c}.multi.txt"))
        per.append(pd.DataFrame({"cohort": c, "default": d, "multi": m}))
    cnt = pd.concat(per).rename_axis("gene_id").reset_index()
    cnt = cnt[cnt["gene_id"].isin(prim.index)]
    cnt["alt_copy"] = cnt["gene_id"].map(prim["alt_copy"])
    cnt["ratio"] = cnt["multi"] / cnt["default"].replace(0, np.nan)

    summary = {"primary_genes": int(len(prim)), "alt_copy_primary_genes": int(prim["alt_copy"].sum()),
               "alt_copy_protein_coding": int((prim["alt_copy"] & (prim["gene_type"] == "protein_coding")).sum()),
               "by_cohort": {}}
    for c, d in cnt.groupby("cohort"):
        e = d[d["default"] >= a.min_count]
        hid = d[(d["default"] < 10) & (d["multi"] >= a.min_count)]
        summary["by_cohort"][c] = {
            grp: {"genes_tested": int(len(x)),
                  "median_ratio": float(x["ratio"].median()),
                  "frac_ratio_ge_1.5": float((x["ratio"] >= 1.5).mean()),
                  "frac_ratio_ge_2": float((x["ratio"] >= 2).mean()),
                  "hidden_genes_default_lt10_multi_ge50": int(hid["alt_copy"].eq(lab).sum())}
            for grp, lab, x in [("alt_copy", True, e[e["alt_copy"]]), ("other", False, e[~e["alt_copy"]])]}

    wide = cnt.pivot_table(index="gene_id", columns="cohort", values="ratio")
    wide["mean_ratio"] = wide.mean(axis=1)
    wide = wide.join(prim[["gene_name", "gene_type", "alt_copy"]])

    deg = pd.read_csv(a.degs)
    deg["canonical_deg"] = (deg["padj"] < 0.05) & (deg["logFC"].abs() > 0.5)
    wide = wide.join(deg.set_index("gene")[["logFC", "padj", "canonical_deg"]], how="left")
    wide["in_deg_table"] = wide["padj"].notna()
    cd = wide[wide["canonical_deg"] == True]  # noqa: E712
    summary["canonical_degs"] = {
        "total": int(deg["canonical_deg"].sum()),
        "alt_copy_among_degs": int(cd["alt_copy"].sum()),
        "degs_mean_ratio_ge_1.5": int((cd["mean_ratio"] >= 1.5).sum()),
        "alt_copy_degs_mean_ratio_ge_1.5": int(((cd["mean_ratio"] >= 1.5) & cd["alt_copy"]).sum()),
    }
    summary["alt_copy_genes_absent_from_deg_table"] = int((~wide.loc[wide["alt_copy"] == True, "in_deg_table"]).sum())  # noqa: E712
    wide.sort_values("mean_ratio", ascending=False).to_csv(out / "gene_multimapper_ratios.tsv", sep="\t")
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
