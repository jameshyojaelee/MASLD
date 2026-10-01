#!/usr/bin/env python3
"""Model A, T4 for tag set v2: allele depths at every v2 tag SNV.

v2 tag SNVs are the union of the v2 `tags/chr*.tags.tsv` positions. Positions shared with
v1 take their rows from the v1 long table (same pileup settings); positions new in v2 take
their rows from the v2-new collect output. Fails if any v2 SNV is in neither source or if a
position appears in both with different REF/ALT. v1 positions dropped from v2 are left out.

Also writes an outcome-free depth summary (total depth n = ref_reads + alt_reads only; no
allele fraction is computed or written) per cohort, for sizing the likelihood computation.
The long table stays in the restricted dir; the depth summary is aggregate.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

DEPTH_CUTS = [1, 8, 20, 50, 100, 200, 500, 1000, 2000, 5000]
QUANTS = [0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags-dir", required=True, help="v2 tags/ directory with chr*.tags.tsv")
    ap.add_argument("--v1-ad", required=True)
    ap.add_argument("--v2new-ad", required=True)
    ap.add_argument("--out", required=True, help="merged long table (restricted dir)")
    ap.add_argument("--depth-summary", required=True, help="aggregate depth table (tsv)")
    a = ap.parse_args()

    tags = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(Path(a.tags_dir).glob("chr*.tags.tsv"))])
    snv = tags.drop_duplicates(["chrom", "pos"])[["chrom", "pos", "ref", "alt"]]
    key = ["chrom", "pos", "ref", "alt"]
    v1 = pd.read_csv(a.v1_ad, sep="\t")
    new = pd.read_csv(a.v2new_ad, sep="\t")

    v1_pos = v1[key].drop_duplicates()
    new_pos = new[key].drop_duplicates()
    both = v1_pos.merge(new_pos, on=["chrom", "pos"], suffixes=("_v1", "_new"))
    clash = both[(both["ref_v1"] != both["ref_new"]) | (both["alt_v1"] != both["alt_new"])]
    if len(clash):
        raise SystemExit(f"{len(clash)} positions differ in REF/ALT between v1 and v2-new tables")

    src = snv.merge(v1_pos.assign(in_v1=True), on=key, how="left").merge(
        new_pos.assign(in_new=True), on=key, how="left")
    src[["in_v1", "in_new"]] = src[["in_v1", "in_new"]].fillna(False).astype(bool)
    # a v2 SNV with no read in any library of any cohort is absent from both tables; count it
    missing = src[~src["in_v1"] & ~src["in_new"]]

    take_v1 = v1.merge(snv, on=key)
    take_new = new.merge(src.loc[~src["in_v1"], key], on=key)
    out = pd.concat([take_v1, take_new], ignore_index=True)
    if out.duplicated(["cohort", "run", "chrom", "pos"]).any():
        raise SystemExit("duplicate library x SNV rows after merge")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False, compression="gzip")

    n = out["ref_reads"] + out["alt_reads"]
    rows = []
    for c, g in n.groupby(out["cohort"]):
        r = {"cohort": c, "libraries": out.loc[g.index, "run"].nunique(), "rows_n_ge_1": int((g >= 1).sum())}
        r.update({f"q{q}": float(np.quantile(g, q)) for q in QUANTS})
        r.update({f"rows_n_ge_{d}": int((g >= d).sum()) for d in DEPTH_CUTS[1:]})
        r["max_n"] = int(g.max())
        rows.append(r)
    pd.DataFrame(rows).to_csv(a.depth_summary, sep="\t", index=False)

    print(f"v2 tag SNVs {len(snv)}; from v1 table {int(src['in_v1'].sum())}; "
          f"from v2-new table only {int((~src['in_v1'] & src['in_new']).sum())}; "
          f"with no reads in any library {len(missing)}")
    print(f"rows written {len(out)} (v1 {len(take_v1)}, new {len(take_new)}); libraries {out['run'].nunique()}")


if __name__ == "__main__":
    main()
