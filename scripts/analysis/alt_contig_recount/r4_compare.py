#!/usr/bin/env python3
"""R4: compare the F_five DEG table rebuilt from original counts (orig/) and from
ALT-copy-corrected counts (corr/). orig/ must reproduce the adopted F_five table;
otherwise the comparison is not interpretable and the script stops.

Canonical DEG = padj < 0.05 and |logFC| > 0.5 (load_figure_data.R::is_canonical_deg).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def deg(d):
    return (d["padj"] < 0.05) & (d["logFC"].abs() > 0.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reference", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    ref = pd.read_csv(a.reference).set_index("gene")
    orig = pd.read_csv(out / "roots/orig/results/integration/deg_results.csv").set_index("gene")
    corr = pd.read_csv(out / "roots/corr/results/integration/deg_results.csv").set_index("gene")
    rep = pd.read_csv(out / "roots/replaced_gene_counts.tsv.gz", sep="\t")
    alt_ids = set(rep["gene_id"])

    shared = ref.index.intersection(orig.index)
    repro = {
        "reference_genes": int(len(ref)), "orig_genes": int(len(orig)),
        "same_gene_set": bool(set(ref.index) == set(orig.index)),
        "max_abs_logFC_diff": float((ref.loc[shared, "logFC"] - orig.loc[shared, "logFC"]).abs().max()),
        "reference_degs": int(deg(ref).sum()), "orig_degs": int(deg(orig).sum()),
    }
    repro["reproduced"] = bool(repro["same_gene_set"] and repro["max_abs_logFC_diff"] < 1e-8
                               and repro["reference_degs"] == repro["orig_degs"])
    summary = {"reproduction": repro}
    if not repro["reproduced"]:
        (out / "comparison.json").write_text(json.dumps(summary, indent=2))
        raise SystemExit("orig/ did not reproduce the adopted F_five table; comparison withheld")

    o_deg, c_deg = set(orig.index[deg(orig)]), set(corr.index[deg(corr)])
    newly_tested = corr.index.difference(orig.index)
    dropped = orig.index.difference(corr.index)
    both = orig.index.intersection(corr.index)
    alt_both = both.intersection(pd.Index(sorted(alt_ids)))
    d_lfc = (corr.loc[alt_both, "logFC"] - orig.loc[alt_both, "logFC"])
    summary.update({
        "genes_tested": {"orig": int(len(orig)), "corr": int(len(corr)),
                         "newly_tested": int(len(newly_tested)), "no_longer_tested": int(len(dropped)),
                         "newly_tested_alt_copy": int(len(newly_tested.intersection(pd.Index(sorted(alt_ids)))))},
        "canonical_degs": {"orig": len(o_deg), "corr": len(c_deg),
                           "orig_up": int((deg(orig) & (orig["logFC"] > 0)).sum()),
                           "corr_up": int((deg(corr) & (corr["logFC"] > 0)).sum()),
                           "gained": len(c_deg - o_deg), "lost": len(o_deg - c_deg),
                           "gained_alt_copy": len((c_deg - o_deg) & alt_ids),
                           "lost_alt_copy": len((o_deg - c_deg) & alt_ids),
                           "gained_non_alt": len((c_deg - o_deg) - alt_ids),
                           "lost_non_alt": len((o_deg - c_deg) - alt_ids)},
        "alt_copy_genes_tested_in_both": int(len(alt_both)),
        "alt_copy_abs_logFC_change_quantiles": {str(q): float(v) for q, v in d_lfc.abs().quantile([0.5, 0.9, 0.99]).items()},
        "non_alt_max_abs_logFC_change": float((corr.loc[both.difference(alt_both), "logFC"] - orig.loc[both.difference(alt_both), "logFC"]).abs().max()),
    })
    cols = ["symbol", "AveExpr", "logFC", "padj"]
    table = orig[cols].add_prefix("orig_").join(corr[cols].add_prefix("corr_"), how="outer")
    table["alt_copy"] = table.index.isin(alt_ids)
    table["orig_deg"] = table.index.isin(o_deg)
    table["corr_deg"] = table.index.isin(c_deg)
    table["symbol"] = table["corr_symbol"].fillna(table["orig_symbol"])
    changed = table[(table["orig_deg"] != table["corr_deg"]) | (table["alt_copy"] & table["orig_logFC"].isna())]
    changed.to_csv(out / "deg_changes.tsv", sep="\t")
    table[table["alt_copy"]].sort_values("corr_AveExpr", ascending=False).to_csv(out / "alt_copy_genes_orig_vs_corr.tsv", sep="\t")
    (out / "comparison.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
