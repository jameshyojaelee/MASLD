"""Fig 2F geometry: map each data circle to (gene row, ancestry column); fit diameter = k*sqrt(value);
compare to recomputed shown values and states. usage: geom_F.py <dump prefix> <tag>"""
import sys, json, numpy as np, pandas as pd
pre, tag = sys.argv[1], sys.argv[2]
cv = pd.read_csv(f"{pre}.curves.tsv", sep="\t"); wd = pd.read_csv(f"{pre}.words.tsv", sep="\t")
cols = wd[wd.text.isin(["EUR", "AFR", "AMR", "EAS", "SAS"]) & (wd.top < 10)].copy(); cols["cx"] = (cols.x0 + cols.x1) / 2
fc = pd.read_csv("out_ABCEF/F_cells.tsv", sep="\t")
genes = set(fc.symbol)
rows = wd[wd.text.isin(genes) & (wd.x1 < 40)].copy(); rows["cy"] = (rows.top + rows.bottom) / 2
legend_top = wd[wd.text.isin(["multi-signal", "PP.H4"])].top.min()
dat = cv[cv.top < legend_top - 2].copy()
dat["cx"] = (dat.x0 + dat.x1) / 2; dat["cy"] = dat.top
dat["col"] = dat.cx.apply(lambda x: cols.iloc[(cols.cx - x).abs().argmin()].text)
dat["row"] = dat.cy.apply(lambda y: rows.iloc[(rows.cy - y).abs().argmin()].text)
dat["drow"] = dat.cy.apply(lambda y: (rows.cy - y).abs().min())
m = dat.merge(fc, left_on=["row", "col"], right_on=["symbol", "anc5"], how="left")
# legend size key
leg = cv[(cv.top >= legend_top - 2) & (cv.fill == "#FFFFFF")].sort_values("w")
legw = leg.w.tolist()
legv = [0.25, 0.5, 0.75, 1.0]
k_leg = float(np.mean(np.array(legw) / np.sqrt(legv))) if len(legw) == 4 else None
# fit to data, excluding tiny marks (minimum-size floor)
ok = m.shown.notna() & (m.w > 1.01)
k_dat = float(np.sum(m.w[ok] * np.sqrt(m.shown[ok])) / np.sum(m.shown[ok]))
m["implied_value_legend_scale"] = (m.w / k_leg) ** 2 if k_leg else np.nan
m["implied_value_data_scale"] = (m.w / k_dat) ** 2
colmap = {"#1565C0": "multi", "#64B5F6": "single_only", "#B0BEC5": "le_0.5"}
m["drawn_state"] = m.fill.map(colmap)
m["state_ok"] = m.drawn_state == m.state
m["abs_err_legend"] = (m.implied_value_legend_scale - m.shown).abs()
m.to_csv(f"out_geom/F_{tag}_matched.tsv", sep="\t", index=False)
row_order = rows.sort_values("cy").text.tolist()
out = {"tag": tag, "n_data_marks": len(dat), "cells_recomputed": len(fc), "unmatched_cells": int(m.shown.isna().sum()),
       "max_row_dist": float(dat.drow.max()), "legend_diams": legw, "k_legend": k_leg, "k_data_fit": k_dat,
       "state_mismatches": m[~m.state_ok][["row", "col", "fill", "state", "s4", "a4"]].to_dict("records"),
       "n_value_off_gt_0.05_legend_scale": int((m.abs_err_legend > 0.05).sum()),
       "max_abs_err_legend_scale": float(m.abs_err_legend.max()),
       "drawn_row_order": row_order, "fill_counts": dat.fill.value_counts().to_dict(),
       "tiny_marks": m[m.w <= 1.01][["row", "col", "shown", "state", "w"]].to_dict("records")}
json.dump(out, open(f"out_geom/F_{tag}.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
print(m[["row", "col", "state", "drawn_state", "shown", "w", "implied_value_legend_scale"]].sort_values("abs_err_legend" if False else "w").tail(12).to_string(index=False) if False else "")
print(m.assign(err=m.abs_err_legend).sort_values("err", ascending=False).head(12)[["row", "col", "state", "drawn_state", "s4", "a4", "shown", "w", "implied_value_legend_scale"]].to_string(index=False))
