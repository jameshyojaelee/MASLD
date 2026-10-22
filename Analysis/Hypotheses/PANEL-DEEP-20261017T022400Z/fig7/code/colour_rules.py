"""Check each panel's colour encoding against the rule its own legend states, using this review's recomputed
intervals (not the drawn ones). Usage: python colour_rules.py <out_dir> <main|calib>"""
import sys
import pandas as pd
OUT, PRE = sys.argv[1], sys.argv[2]
t = pd.read_csv(f"{OUT}/geometry_{PRE}.tsv", sep="\t")
MAG, BLU, GRY = "#C9265E", "#1565C0", "#9E9E9E"
rows = []
# 7B: above zero magenta, below zero blue, crosses zero gray (legend), applied to the drawn 95% interval
b = t[t.panel == "fig7B_gain_over_task_native_baseline"]
for i, r in b[b.mark == "point"].iterrows():
    sub = b[(b.index > i) & (b.index <= i + 2)]
    lo = sub[sub.mark == "interval_lo"].recomputed; hi = sub[sub.mark == "interval_hi"].recomputed
    if len(lo) and len(hi):
        lo, hi = float(lo.iloc[0]), float(hi.iloc[0])
        exp = MAG if lo > 0 else BLU if hi < 0 else GRY
        icol = sub[sub.mark == "interval_lo"].colour.iloc[0]
    else:
        exp = GRY; icol = ""
    pcol = r.colour.split("/")[1]
    rows.append({"panel": "7B", "row": r.row, "rule": "above/below/crosses zero (95% interval)", "expected": exp, "point": pcol, "interval": icol,
                 "ok": pcol == exp and (icol in ("", exp))})
# 7F gains: above zero magenta fill, crosses zero gray, no clean study hollow
f = t[t.panel == "fig7F_liver_encoder_against_the_field"]
for i, r in f[f.mark == "point"].iterrows():
    sub = f[(f.index > i) & (f.index <= i + 2)]
    lo = float(sub[sub.mark == "interval_lo"].recomputed.iloc[0]); hi = float(sub[sub.mark == "interval_hi"].recomputed.iloc[0])
    pfill = r.colour.split("/")[1]
    exp = "none" if r.row.startswith("geneformer") else (MAG if lo > 0 else BLU if hi < 0 else GRY)
    icol = sub[sub.mark == "interval_lo"].colour.iloc[0]
    exp_i = GRY if exp == "none" else exp
    rows.append({"panel": "7F", "row": r.row, "rule": "above zero / crosses zero / no clean study", "expected": exp, "point": pfill, "interval": icol,
                 "ok": pfill == exp and icol == exp_i})
# 7E: colour by verdict group
e = t[(t.panel == "fig7E_verdict_vocabulary_and_refusal_thresholds") & (t.mark == "point")]
for _, r in e.iterrows():
    v = r.note.replace("verdict drawn: ", "")
    exp = {"supported": MAG, "untestable": GRY, "tested negative": BLU}[v]
    rows.append({"panel": "7E", "row": r.row, "rule": "colour = verdict group", "expected": exp, "point": r.colour, "interval": "", "ok": r.colour == exp})
# 7C: exposure resolved magenta; training-seen orange
c = t[(t.panel == "fig7C_score_beside_training_exposure") & (t.mark == "point")]
for _, r in c.iterrows():
    exp = "#F4A674" if r.row == "TranscriptFormer" else MAG
    rows.append({"panel": "7C", "row": r.row, "rule": "exposure resolved / includes training-seen", "expected": exp, "point": r.colour.split("/")[1], "interval": "", "ok": r.colour.split("/")[1] == exp})
o = pd.DataFrame(rows)
o.to_csv(f"{OUT}/colour_rules_{PRE}.tsv", sep="\t", index=False)
print(PRE, "checked", len(o), "violations", int((~o.ok).sum()))
print(o[~o.ok].to_string())
