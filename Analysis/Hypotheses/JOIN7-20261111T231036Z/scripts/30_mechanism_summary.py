#!/usr/bin/env python3
"""JOIN7 step 3 (after job 23226934): one row per mechanism from mechanism_endpoint_table.tsv.
Usage: 30_mechanism_summary.py <out_dir>   Writes mechanism_summary.tsv and mechanism_summary.md."""
import sys
import pandas as pd

out = sys.argv[1]
d = pd.read_csv(f"{out}/mechanism_endpoint_table.tsv", sep="\t", dtype=str)
SUP = {"gene_level_support", "directional_agree", "directional_disagree", "mixed_direction"}
EP = ["pdff", "alt", "mash_resolution", "fibrosis_improvement"]
SHORT = {"improved": "improved", "improved_unadjusted": "improved (unadj.)", "not_improved": "not improved",
         "no_test_reported": "not reported"}


def sup(row, t):
    c = row[f"class_{t}"]
    if c not in SUP:
        return c.replace("_", " ")
    genes = [x.split(":")[0] for x in row[f"genes_{t}"].split(";") if x.split(":")[1] in SUP]
    return "support (" + ",".join(genes) + ")"


rows = []
for mech, g in d.groupby("mechanism", sort=False):
    r0 = g.iloc[0]
    rec = dict(mechanism=mech, targets=r0.targets_used if isinstance(r0.targets_used, str) else "none",
               action=r0.action_direction, liver_fat=sup(r0, "pdff"), cirrhosis=sup(r0, "cirr"),
               alt_ukb=sup(r0, "alt_ukb"), alt_mvp=sup(r0, "alt_mvp"), direction_call="none mapped",
               not_verified_rows=int((g.not_verified == "True").sum()))
    for e in EP:
        x = g[g.endpoint == e]
        rec[e] = "" if x.empty else f"{SHORT[x.result_vs_control.iloc[0]]}; {x.drug.iloc[0]} {x.nct_id.iloc[0]} ph{x.phase.iloc[0]}"
    rows.append(rec)
m = pd.DataFrame(rows)
m.to_csv(f"{out}/mechanism_summary.tsv", sep="\t", index=False)
cols = ["mechanism", "targets", "liver_fat", "cirrhosis", "alt_ukb"] + EP + ["not_verified_rows"]
with open(f"{out}/mechanism_summary.md", "w") as f:
    f.write("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n")
    for _, r in m.iterrows():
        f.write("| " + " | ".join(str(r[c]) if str(r[c]) else "-" for c in cols) + " |\n")
print(len(m), "mechanisms")
