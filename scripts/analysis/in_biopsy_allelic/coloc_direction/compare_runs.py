#!/usr/bin/env python3
"""B-COLOC: which gene x release x trait calls changed between two direction runs, and the stated reasons.

Usage: compare_runs.py <old run dir> <new run dir> <out tsv>
Joins direction_by_trait.tsv of both runs on release x ensembl x trait and writes one row per
row of either run; prints per-release counts (JSON). The cause of a change is read from the
not_directional_reason of the side that is not directional (new side for directional ->
not_directional, old side for not_directional -> directional).
"""
import json
import sys
from pathlib import Path

import pandas as pd

KEY = ["release", "ensembl", "trait"]
COLS = ["gene", "direction_sign", "not_directional_reason", "reported_study", "hit1", "hit2"]


def load(d):
    t = pd.read_csv(Path(d) / "direction_by_trait.tsv", sep="\t", dtype=str, keep_default_na=False)
    t = t.rename(columns={"driving_study": "reported_study"})      # v1 (full-20260926T171722Z) name
    t["direction_sign"] = pd.to_numeric(t.direction_sign).astype(int)
    return t[KEY + COLS]


def norm(s):
    return s.str.replace(r" in (EUR|EAS|AFR|AMR|SAS) ", " in <pop> ", regex=True) \
            .str.replace(r"\(([^)]*)\)$", "", regex=True).str.strip()


def main():
    old_dir, new_dir, out = sys.argv[1:4]
    m = load(old_dir).merge(load(new_dir), on=KEY, how="outer", suffixes=("_old", "_new"), indicator=True)
    so, sn = m.direction_sign_old, m.direction_sign_new
    m["change"] = "same"
    m.loc[(so != 0) & (sn != 0) & (so != sn), "change"] = "sign flip"
    m.loc[(so != 0) & (sn == 0), "change"] = "directional -> not_directional"
    m.loc[(so == 0) & (sn != 0), "change"] = "not_directional -> directional"
    m.loc[m._merge != "both", "change"] = "row in one run only"
    m["cause"] = ""
    d2n = m.change == "directional -> not_directional"
    n2d = m.change == "not_directional -> directional"
    m.loc[d2n, "cause"] = "new: " + norm(m.loc[d2n, "not_directional_reason_new"])
    m.loc[n2d, "cause"] = "old: " + norm(m.loc[n2d, "not_directional_reason_old"])
    m = m.drop(columns="_merge").sort_values(["release", "change", "gene_new", "trait"])
    m.to_csv(out, sep="\t", index=False)
    summary = {"old": str(old_dir), "new": str(new_dir), "rows_matched": int((m.change != "row in one run only").sum()),
               "by_release": {}}
    for rel, g in m.groupby("release"):
        summary["by_release"][rel] = {
            "changes": g.change.value_counts().to_dict(),
            "sign_flips": [f"{r.gene_new} {r.trait} {int(r.direction_sign_old):+d} -> {int(r.direction_sign_new):+d}"
                           for r in g[g.change == "sign flip"].itertuples()],
            "directional -> not_directional, new reason": g.loc[g.change == "directional -> not_directional", "cause"]
            .value_counts().to_dict(),
            "not_directional -> directional, old reason": g.loc[g.change == "not_directional -> directional", "cause"]
            .value_counts().to_dict()}
    print(json.dumps(summary, indent=2))
    Path(out).with_suffix(".json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
